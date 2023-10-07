#!/opt/conda/bin/python
import sys
import timeit
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import pandas as pd
from sklearn.metrics import roc_auc_score
import preprocess as pp
import pickle
import random
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import median_absolute_error,r2_score, mean_absolute_error,mean_squared_error

class MolecularGraphNeuralNetwork(nn.Module):
    def __init__(self, N, dim, layer_hidden, layer_output):
        super(MolecularGraphNeuralNetwork, self).__init__()
        self.embed_fingerprint = nn.Embedding(N, dim)
        self.W_fingerprint = nn.ModuleList([nn.Linear(dim, dim)
                                            for _ in range(layer_hidden)])
        self.W_output = nn.ModuleList([nn.Linear(dim, dim)
                                       for _ in range(layer_output)])
        self.W_property = nn.Linear(dim, 1)

    def pad(self, matrices, pad_value):
        """Pad the list of matrices
        with a pad_value (e.g., 0) for batch processing.
        For example, given a list of matrices [A, B, C],
        we obtain a new matrix [A00, 0B0, 00C],
        where 0 is the zero (i.e., pad value) matrix.
        """
        shapes = [m.shape for m in matrices]
        M, N = sum([s[0] for s in shapes]), sum([s[1] for s in shapes])
        zeros = torch.FloatTensor(np.zeros((M, N))).to(device)
        pad_matrices = pad_value + zeros
        i, j = 0, 0
        for k, matrix in enumerate(matrices):
            m, n = shapes[k]
            pad_matrices[i:i+m, j:j+n] = matrix
            i += m
            j += n
        return pad_matrices

    def update(self, matrix, vectors, layer):
        hidden_vectors = torch.relu(self.W_fingerprint[layer](vectors))
        return hidden_vectors + torch.matmul(matrix, hidden_vectors)

    def sum(self, vectors, axis):
        sum_vectors = [torch.sum(v, 0) for v in torch.split(vectors, axis)]
        return torch.stack(sum_vectors)

    def mean(self, vectors, axis):
        mean_vectors = [torch.mean(v, 0) for v in torch.split(vectors, axis)]
        return torch.stack(mean_vectors)

    def gnn(self, inputs):

        """Cat or pad each input data for batch processing."""
        Smiles,fingerprints, adjacencies, molecular_sizes = inputs
        fingerprints = torch.cat(fingerprints)
        adjacencies = self.pad(adjacencies, 0)

        """GNN layer (update the fingerprint vectors)."""
        fingerprint_vectors = self.embed_fingerprint(fingerprints)
        for l in range(layer_hidden):
            hs = self.update(adjacencies, fingerprint_vectors, l)
            fingerprint_vectors = F.normalize(hs, 2, 1)  # normalize.

        """Molecular vector by sum or mean of the fingerprint vectors."""
        molecular_vectors = self.sum(fingerprint_vectors, molecular_sizes)
        return Smiles,molecular_vectors

    def mlp(self, vectors):
        """ regressor based on multilayer perceptron."""
        for l in range(layer_output):
            vectors = torch.relu(self.W_output[l](vectors))
        outputs = self.W_property(vectors)
        return outputs
    def forward_regressor(self, data_batch, train):

        inputs = data_batch[:-1]
        correct_values = torch.cat(data_batch[-1])

        if train:
            Smiles,molecular_vectors = self.gnn(inputs)
            predicted_values = self.mlp(molecular_vectors)
            loss = F.mse_loss(predicted_values, correct_values)
            return loss
        else:
            with torch.no_grad():
                Smiles,molecular_vectors = self.gnn(inputs)
                predicted_values = self.mlp(molecular_vectors)
            predicted_values = predicted_values.to('cpu').data.numpy()
            correct_values = correct_values.to('cpu').data.numpy()
            predicted_values = np.concatenate(predicted_values)
            correct_values = np.concatenate(correct_values)
            return Smiles,predicted_values, correct_values
    def forward_predict(self, data_batch):

            inputs = data_batch
            try:
                Smiles,molecular_vectors = self.gnn(inputs)
            except RuntimeError:
                print(data_batch)
            predicted_values = self.mlp(molecular_vectors)
            predicted_values = predicted_values.to('cpu').data.numpy()
            predicted_values = np.concatenate(predicted_values)
            
            return Smiles,predicted_values


class Trainer(object):
    def __init__(self, model):
        self.model = model
        self.optimizer = optim.Adam(self.model.parameters(), lr=lr)

    def train(self, dataset):
        np.random.shuffle(dataset)
        N = len(dataset)
        loss_total = 0
        for i in range(0, N, batch_train):
            data_batch = list(zip(*dataset[i:i+batch_train]))
            loss = self.model.forward_regressor(data_batch, train=True)
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()
            loss_total += loss.item()
        return loss_total

class Trainer_tf(object):
    def __init__(self, model):
        ##转移学习的optimizer为何不一样？
        self.model = model
        self.optimizer = optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=lr)

    def train(self, dataset):
        np.random.shuffle(dataset)
        N = len(dataset)
        loss_total = 0
        for i in range(0, N, batch_train):
            data_batch = list(zip(*dataset[i:i+batch_train]))
            loss = self.model.forward_regressor(data_batch, train=True)
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()
            loss_total += loss.item()
        return loss_total


class Tester(object):
    def __init__(self, model):
        self.model = model
    def test_regressor(self, dataset):
        N = len(dataset)
        SMILES, Ts, Ys = '', [], []
        SAE = 0  # sum absolute error.
        for i in range(0, N, batch_test):
            data_batch = list(zip(*dataset[i:i+batch_test]))
            (Smiles,  predicted_values,correct_values) = self.model.forward_regressor(
                                               data_batch, train=False)
            SMILES += ' '.join(Smiles) + ' '
            Ts.append(correct_values)
            Ys.append(predicted_values)
            
            SAE += sum(np.abs(predicted_values-correct_values))
        SMILES = SMILES.strip().split()
        T, Y = map(str, np.concatenate(Ts)), map(str, np.concatenate(Ys))
        predictions = '\n'.join(['\t'.join(x) for x in zip(SMILES, T, Y)])
        MAEs = SAE / N  # mean absolute error.
        return MAEs,predictions

    def test_predict(self, dataset):
        N = len(dataset)
        SMILES, Ts, Ys = '', [], []
        SAE = 0  # sum absolute error.
        for i in range(0, N, batch_test):
            data_batch = list(zip(*dataset[i:i+batch_test]))
            (Smiles,  predicted_values) = self.model.forward_predict(data_batch)
            SMILES += ' '.join(Smiles) + ' '
            Ys.append(predicted_values)
        SMILES = SMILES.strip().split()
        Y = map(str, np.concatenate(Ys))
        predictions = '\n'.join(['\t'.join(x) for x in zip(SMILES, Y)])
        return predictions

    def save_MAEs(self, MAEs, filename):
        with open(filename, 'a') as f:
            f.write(MAEs + '\n')
    def save_predictions(self, predictions, filename, ct = False):
        with open(filename, 'w') as f:
            if ct:
                f.write('Smiles\tCorrect\tPredict\n')
            else:
                f.write('Smiles\tPredict\n')
            f.write(predictions + '\n')

    def save_model(self, model, filename):
        torch.save(model.state_dict(), filename)

def split_dataset(dataset, ratio):
#   """Shuffle and split a dataset."""
    np.random.seed(1234)  # fix the seed for shuffle.
    np.random.shuffle(dataset)
    n = int(ratio * len(dataset))
    return dataset[:n], dataset[n:]

def dump_dictionary(dictionary, filename):
        with open(filename, 'wb') as f:
            pickle.dump(dict(dictionary), f)

def plot_loss(file_MAEs,test = True):    
    # def rmse(y_true, y_pred):
    #     return np.sqrt(mean_squared_error(y_true, y_pred))
    loss = pd.read_table(file_MAEs)
    plt.plot(loss['MAE_train'], color='r',label='MSE of train set')
    plt.plot(loss['MAE_dev'], color='b',label='MSE of validation set')
    if test:
        plt.plot(loss['MAE_test'], color='y',label='MSE of test set')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend()
    plt.savefig(outp+'loss.png',dpi=300)
    plt.close()

def plot_cp(file_test_result,minutes = True):

    res_tf = pd.read_table(file_test_result)
    r2 = r2_score(res_tf['Correct'], res_tf['Predict'])
    mae = mean_absolute_error(res_tf['Correct'], res_tf['Predict'])
    medae = median_absolute_error(res_tf['Correct'], res_tf['Predict'])
    rmae = np.mean(np.abs(res_tf['Correct'] - res_tf['Predict']) / res_tf['Correct']) * 100
    median_re = np.median(np.abs(res_tf['Correct'] - res_tf['Predict']) / res_tf['Correct'])
    mean_re=np.mean(np.abs(res_tf['Correct'] - res_tf['Predict']) / res_tf['Correct'])
    rtlim = 23 if minutes else 1400
    plt.plot(res_tf['Correct'], res_tf['Predict'], '.', color = 'blue')
    plt.plot([0,rtlim], [0,rtlim], color ='red')
    plt.ylabel('Predicted RT')
    plt.xlabel('Experimental RT')        
    plt.text(0,0.95*rtlim, 'R2='+str(round(r2,4)), fontsize=12)
    plt.text(0.36*rtlim,0.95*rtlim,'MAE='+str(round(mae,4)),fontsize=12)
    plt.text(0, 0.85*rtlim, 'MedAE='+str(round(medae,4)), fontsize=12)
    plt.text(0.36*rtlim, 0.85*rtlim, 'MRE='+str(round(mean_re,4)), fontsize=12)
    plt.text(0, 0.72*rtlim, 'MedRE='+str(round(median_re,4)), fontsize=12)
    plt.savefig( outp+'c-p.png',dpi=300)
    plt.close()


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = None)
    parser.add_argument('-n','--name',help = "name of the orignial model_training dataset ",type = str, default = "SMRT")
    parser.add_argument('-i','--inputname',help = "name of the using dataset for model_training",type = str, default = "BP")
    parser.add_argument('-t','--train',help = "input datasets for trainning",type = str, default = "BP_train_set.csv")
    parser.add_argument('-p','--pred',help = "input datasets for testing",type = str, default = "HMDB_pred_set.csv")
    parser.add_argument('-d','--dir', help = "path for input & output datasets (with slash)",type = str, default = "test/test2/")
    parser.add_argument('-m','--model', help = "input model name",type = str, default = "SMRT_model.h5")
    parser.add_argument('-o','--out', help = "output path",type = str, default = "result/test2/")
    # args = parser.parse_args(["--train","trainset.csv","--pred","predset.csv","--dir","test/test2/data/", "--out", "test/test2/result/"])
    
    args = parser.parse_args()
    
    radius=1
    dim=48
    layer_hidden=6
    layer_output=6
    batch_train=32
    batch_test=32
    lr=2e-4
    lr_decay=0.99
    decay_interval=10
    iteration_tf=10
    N=5000

    #训练数据集
    path= args.dir
    outp = args.out
    data_name=args.name
    idata_name = args.inputname
    train_file = args.train
    pred_file = args.pred
    model_name = path+args.model

    if torch.cuda.is_available():
        device = torch.device('cuda')
        print('The code uses a GPU!')
    else:
        device = torch.device('cpu')
        print('The code uses a CPU...')
   
    print('The code uses a CPU!')

    data = pd.read_csv(path+train_file,sep = "\t")
    
    x=list(data['smiles'])
    y=list(data['RT'])
    index = [i for i in range(len(x))]
    random.shuffle(index)
    X,Y=[],[]
    k=10
    for i in range(len(index)):
        X.append(x[index[i]])
        Y.append(y[index[i]])

    def get_k_fold_data(k, i, X, y): 
        #
        assert k > 1
        fold_size = len(X) // k  
        X_train, Y_train = [], []
        for j in range(k):
            if j==k-1:
                idx = [j * fold_size, len(X)]
            else:
                idx = [j * fold_size, (j + 1) * fold_size] 
        
            X_part, Y_part = X[idx[0]: idx[1]], Y[idx[0]:idx[1]]
            if j == i:
                X_test, Y_test = X_part, Y_part
            elif len(X_train)  == 0:
                X_train, Y_train = X_part, Y_part
            else:
                X_train=X_train+X_part
                Y_train=Y_train+Y_part
        
        dataset_tf_test = pp.create_dataset_kfold(X_test,Y_test,path,outp,data_name,idata_name)

        dataset_tf_train = pp.create_dataset_kfold(X_train,Y_train,path,outp,data_name,idata_name)
        dataset_tf_train, dataset_tf_dev = split_dataset(dataset_tf_train, 0.9)


        return dataset_tf_train, dataset_tf_dev,dataset_tf_test 
    dataset_tf_train, dataset_tf_dev,dataset_tf_test  = get_k_fold_data(k, 0, X, Y) 
    print('-'*100)
    print('# of training data samples:', len(dataset_tf_train))
    print('# of development data samples:', len(dataset_tf_dev))
    print('# of test data samples:', len(dataset_tf_test))

    if pred_file:
        pred_data = pd.read_csv(path + pred_file,sep = "\t", usecols = ["smiles"])
        pred_data = [x for x in list(pred_data.smiles) if "." not in x]
        dataset_tf_pred = pp.transferlearning_dataset_predict(pred_data,path,data_name)
        print('# of prediction data samples:', len(dataset_tf_pred))
    
    print('The preprocess has finished!')
    print('-'*100)
    print('Creating a model.')
    torch.manual_seed(1234)
    model= MolecularGraphNeuralNetwork(N, dim, layer_hidden, layer_output).to(device)

    ##加载已有模型
    model.load_state_dict(torch.load(model_name, map_location=torch.device('cpu')))
    for para in model.W_fingerprint.parameters():
        para.requires_grad = False   
    print(model)

    trainer = Trainer_tf(model)
    tester = Tester(model)
    print('# of model parameters:', sum([np.prod(p.size()) for p in model.parameters()]))
    print('-'*100)

    file_MAEs = outp +'MAEs'+'.txt'
    file_test_result  =  outp + 'test_prediction'+ '.txt'
    file_predictions =  outp + 'prediction' +'.txt'
    file_model =   outp + idata_name +'_model'+'.h5'
    result_tf = 'Epoch\tTime(sec)\tLoss_train\tMAE_train\tMAE_dev\tMAE_test'

    with open(file_MAEs, 'w') as f:
        f.write(result_tf + '\n')
    print('Start training.')
    print('The result is saved in the output directory every epoch!')

    np.random.seed(1234)  
    start = timeit.default_timer()
    for epoch in range(iteration_tf):
        epoch += 1
        if epoch % decay_interval == 0:
            trainer.optimizer.param_groups[0]['lr'] *= lr_decay
        model.train()
        loss_train = trainer.train(dataset_tf_train)
        MAE_tf_best=9999999
        model.eval()
        MAE_tf_train,predictions_train_tf = tester.test_regressor(dataset_tf_train)
        MAE_tf_dev = tester.test_regressor(dataset_tf_dev)[0]
        MAE_tf_test = tester.test_regressor(dataset_tf_test)[0]
        predictions_test_tf  = tester.test_regressor(dataset_tf_test)[1]
        if pred_file:
            predictions_pred_tf = tester.test_predict(dataset_tf_pred)
            tester.save_predictions(predictions_pred_tf , file_predictions, ct = False)

        tester.save_predictions(predictions_test_tf , file_test_result, ct = True)
        time = timeit.default_timer() - start
        if epoch == 1:
            minutes = time * iteration_tf / 60
            hours = int(minutes / 60)
            minutes = int(minutes - 60 * hours)
            print('The training will finish in about',
                  hours, 'hours', minutes, 'minutes.')
            print('-'*100)
            print(result_tf)

        results_tf = '\t'.join(map(str, [epoch, time, loss_train,MAE_tf_train,MAE_tf_dev, MAE_tf_test]))
        tester.save_MAEs(results_tf, file_MAEs)

        if MAE_tf_dev <= MAE_tf_best:
            MAE_tf_best = MAE_tf_dev
            tester.save_model(model, file_model)
        print(results_tf)
    

    plot_loss(file_MAEs,test = True)        
    plot_cp(file_test_result)
