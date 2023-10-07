#!/opt/conda/bin/python
import email
import time
import os
import re
import datetime
import base64
import pytz
import click
import pickle
import pathlib
from imaplib import IMAP4_SSL
from bs4 import BeautifulSoup
from typing import Union
import pandas as pd
import multiprocessing as mp


@click.command()
@click.option("-u", "--user", type=str, help="like love@qq.com")
@click.option("-p", "--password", type=str, help="password for user")
@click.option("-m", "--mailbox", type=str, help='mailbox such as : "Sent Messages" or INBOX')
@click.option("-n", "--name", type=str, help="such as: 数据")
@click.option("-s", "--since", type=str, help="date like: 30-Nov-2022")
@click.option("-l", "--load", type=bool,default = True, help="if load data")
def main_task(user, password, mailbox, name, since, load=True):
    """
    邮件爬虫主流程
    tool_mail_getinfo --user "zhaobx@lumingbio.com" --password PsDsKhBq47KWAxgm --mailbox '"Sent Messages"' --name "原始数据" --since "01-Jun-2023"
    tool_mail_getinfo --user "zhaobx@lumingbio.com" --password PsDsKhBq47KWAxgm --mailbox '"Sent Messages"' --name "报告" --since "01-Jun-2023"

    tool_mail_getinfo --user "zhouy@lumingbio.com" --password ESEhK2QJk7P42jc5 --mailbox '"Sent Messages"' --name "原始数据" --since "01-Jun-2023"
    tool_mail_getinfo --user "zhouy@lumingbio.com" --password ESEhK2QJk7P42jc5 --mailbox '"Sent Messages"' --name "报告" --since "01-Jun-2023"

    tool_mail_getinfo --user "wanghy@lumingbio.com" --password Kzjy7bmq8eUqxhhp --mailbox '"Sent Messages"' --name "原始数据" --since "01-Jun-2023"
    tool_mail_getinfo --user "wanghy@lumingbio.com" --password Kzjy7bmq8eUqxhhp --mailbox '"Sent Messages"' --name "报告" --since "01-Jun-2023"
    """
    # con1 = [{'name': '数据', 'SINCE': '30-Nov-2022', 'SUBJECT': '数据'}]
    con = [{'name': name, 'SINCE': since, 'SUBJECT': name}]
    login = {'user': user, 'password': password, 'mailbox': mailbox}
    name = re.sub('@.*', '_', user) + mailbox.replace(' ', '').replace('"', '').replace("'", '') + '_'
    p_num = 8
    # 登录邮箱
    email_server = Email(user=user, password=password, mailbox=mailbox)
    email_server.link_server()

    header_log = name + 'header.log'
    body_log = name + 'body.log'
    fail_log_path = name + 'fail_log.xlsx'
    header_tmp_dir = name + 'header_tmp'
    body_tmp_dir = name + 'body_tmp'
    ## 获取邮件头信息
    header_list = []
    # try load header
    if load:
        try:
            for v in range(len(con)):
                header_path = name + str(con[v].get('name')) + '_header.xlsx'
                header_list.append(pd.read_excel(header_path, index_col=None))
                print('load ' + header_path)
        except:
            print('load fail')
            header_list = []
    else:
        if os.path.exists(header_log):
            os.remove(header_log)
        if os.path.exists(body_log):
            os.remove(body_log)
    if os.path.exists(header_log):
        log = read_log(header_log)
        uids = list(log.keys())
    else:
        # get all uids by try search
        ### 将需要搜索的邮件主题进行编码,拼接搜索条件
        search_query = 'OR'
        date = remove_timezone(datetime.datetime.fromtimestamp(time.time()))
        for v in con:
            search_query = search_query + ' SUBJECT "{}"'.format(v.get('SUBJECT').encode('unicode_escape'))
            date_tmp = structure_date(v.get('SINCE'))
            if date_tmp < date:
                date = date_tmp
        ### 使用 uid() 方法进行搜索
        typ, msg_nums = email_server.email_server.uid('SEARCH', 'SINCE "{}"'.format(date.strftime("%d-%b-%Y")),
                                                      search_query)
        uids = msg_nums[0].decode().split()

    if header_list == [] or len(uids) > 0:
        # 抓取邮件头
        arg = '(BODY[HEADER.FIELDS (Subject DATE From To)])'
        if not os.path.exists(header_tmp_dir):
            os.mkdir(header_tmp_dir)
        uids_tmp = check_log(log=uids, res_dir=header_tmp_dir)
        if len(uids_tmp) > 0:
            fetch_mail(uids=uids_tmp, login=login, arg=arg, log_path=header_log, res_dir=header_tmp_dir, p_num=p_num)
        # 解析邮件头
        header_list, fail_log = parse_header(con=con, path=header_tmp_dir, log=uids)
        for v in range(len(con)):
            header_path = name + str(con[v].get('name')) + '_header.xlsx'
            header_list[v].to_excel(header_path, index=None)
        if os.path.exists(fail_log_path):
            os.remove(fail_log_path)
        if fail_log is not None:
            fail_log.to_excel(fail_log_path, index=None)
        if os.path.exists(body_log):
            os.remove(body_log)
    # 抓取邮件主体
    for v in range(len(con)):
        uids += list(header_list[v]['uid'])
    uids = set(uids)
    arg = '(RFC822)'
    if not os.path.exists(body_tmp_dir):
        os.mkdir(body_tmp_dir)
    uids_tmp = check_log(log=list(uids), res_dir=body_tmp_dir)
    if len(uids_tmp) > 0:
        fetch_mail(uids=uids_tmp, login=login, arg=arg, log_path=body_log, res_dir=body_tmp_dir, p_num=p_num)

    # 解析邮件主体 path:error
    res = parse_body(path=body_tmp_dir, log=list(uids))
    res_list = []
    for v in range(len(con)):
        res_path = name + str(con[v].get('name')) + '_total.xlsx'
        res_list.append(header_list[v].merge(res, on='uid', how='left'))
        res_list[v].to_excel(res_path, index=None)


def check_log(log: Union[str, list, pathlib.Path], res_dir=None):
    if res_dir is not None:
        res_list = os.listdir(res_dir)
    else:
        res_list = os.listdir()
    res_list = pd.DataFrame({'id': res_list})
    res_list = res_list.applymap(lambda x: re.sub('_.*', '', str(x)))
    if isinstance(log, list):
        log = pd.DataFrame({'id': log})
    else:
        log = pd.read_csv(log, sep='\t')
    log = log.applymap(lambda x: str(x))
    con = log.iloc[:, 0].isin(res_list["id"]).map(lambda x: not x)
    log = log.loc[con, :]
    uids = list(log.iloc[:, 0])
    return uids


def read_log(log_path):
    res = {}
    log = pd.read_csv(log_path, sep='\t')
    log = log.apply(lambda x: res.update({str(x[0]): str(x[1])}), axis=1)
    return res


def save_data(data, name=None):
    if name is None:
        name = 'data.pickle'
    with open(name, 'wb') as f:
        pickle.dump(data, f)


def read_data(name=None):
    if name is None:
        name = 'data.pickle'
    with open(name, 'rb') as f:
        data = pickle.load(f)
    return data


def mail_dic2pd(res: dict):
    res_list = []
    for uid in res:
        msg_data = res[uid]
        raw_email = get_msg_data(msg_data, name=uid)
        res_list.append({'uid': uid, 'raw': raw_email})
    return pd.DataFrame(res_list)


def get_msg_data(msg_data, name=None):
    if msg_data is None:
        print(str(name) + ' None')
        res = None
    elif isinstance(msg_data, bytes):
        res = msg_data
    elif msg_data[0] is None:
        print(str(name) + ' None')
        res = None
    else:
        res = msg_data[0][1]
    return res


def parse_body(res: dict = None, path=None, log: Union[dict, list] = None):
    res_list = []
    if path is not None:
        path = pathlib.Path(path)
    if isinstance(log, list):
        log_tmp = {}
        for v in log:
            log_tmp.update({v: 0})
        log = log_tmp
    if res is None:
        res = log
    for uid in res:
        if log is not None:
            name = '_'.join([str(uid), str(log[uid])])
            if path is not None:
                name = path / name
            msg_data = read_bytes(path=name)
        else:
            msg_data = res[uid]
        raw_email = get_msg_data(msg_data, name=uid)
        if raw_email is None:
            continue
        soup_list, file_list = message_parsing(raw_email)
        print('deal done')
        res_list.append({'uid': uid, 'text': '\n'.join(soup_list), 'files': ','.join(file_list)})
    data = pd.DataFrame(res_list)
    ILLEGAL_CHARACTERS_RE = re.compile(r'[\000-\010]|[\013-\014]|[\016-\037]')
    data = data.applymap(lambda x: ILLEGAL_CHARACTERS_RE.sub(r'', x) if isinstance(x, str) else x)
    return data


def structure_date(raw_date: str):
    try:
        if isinstance(raw_date, str):
            structure_date = time.strptime(raw_date, "%d-%b-%Y")
            date = datetime.datetime.fromtimestamp(time.mktime(structure_date))
        elif isinstance(raw_date, datetime.datetime):
            date = raw_date
        else:
            raise TypeError('raw_date must be datetime.datetime or %d-%b-%Y')
    except ValueError as e:
        raise ValueError('Please enter the date in format of \'01-Jan-2022\'.', e)
    date = remove_timezone(date)
    return date


def parse_header(con: list, res: dict = None, path=None, log: Union[dict, list] = None):
    """解析邮件头"""
    res_list = []
    fail_log = []
    if path is not None:
        path = pathlib.Path(path)
    if isinstance(log, list):
        log_tmp = {}
        for v in log:
            log_tmp.update({v: 0})
        log = log_tmp
    if res is None:
        res = log
    for v in range(len(con)):
        res_list.append([])
        date = structure_date(con[v].get('SINCE'))
        con[v]['SINCE'] = date
    for uid in res:
        if log is not None:
            name = '_'.join([str(uid), str(log[uid])])
            if path is not None:
                name = path / name
            msg_data = read_bytes(path=name)
        else:
            msg_data = res[uid]
        raw_email = get_msg_data(msg_data, name=uid)
        if raw_email is None:
            fail_log.append(uid)
            continue
        # 将bytes类型的msg_data[0][1]转换为email.message.Message对象，不会出现解码错误。
        msg_header = email.message_from_bytes(raw_email)
        if msg_header['Date'] is None:
            # 没有时间，判断为没有内容
            print(str(uid) + ' None date')
            fail_log.append(uid)
            continue
        # Message对象转换为字符串时可能会出现解码错误。
        # 时间一般没有解码问题，直接读
        Date = email.utils.parsedate_to_datetime(msg_header['Date'])
        Date = remove_timezone(Date)
        # email.header.decode_header 尝试识别字符集类型，可能会出现识别错误
        Subject = decode_header(msg_header, con='Subject', join=',', pattern=None)
        From = decode_header(msg_header, con='From', join=',', pattern='com')
        To = decode_header(msg_header, con='To', join=',', pattern='com')
        for v in range(len(con)):
            if Date > con[v]['SINCE']:
                if re.search(con[v]['SUBJECT'], Subject):
                    print(str(uid) + ' add ' + str(v))
                    res_list[v].append(
                        {'uid': uid, 'Date': Date, 'Subject': Subject, 'From': From, 'To': To, 'text': None,
                         'files': None})
                else:
                    print(str(uid) + ' skip' + str(v))
            else:
                print(str(uid) + ' skip' + str(v))
    data_list = []
    for v in range(len(con)):
        data = pd.DataFrame(res_list[v])
        ILLEGAL_CHARACTERS_RE = re.compile(r'[\000-\010]|[\013-\014]|[\016-\037]')
        data = data.applymap(lambda x: ILLEGAL_CHARACTERS_RE.sub(r'', x) if isinstance(x, str) else x)
        data_list.append(data)
    if len(fail_log) > 0:
        fail_log = pd.DataFrame({'uid': fail_log})
    else:
        fail_log = None
    return data_list, fail_log


def fetch_mail(uids: list, login: dict, arg: str = '(RFC822)', log_path='task_log.log', res_dir='tmp', p_num: int = 4):
    """fetch_mail"""

    def write_tmp(result_tmp: list, res_dir):
        for task in result_tmp:
            res_dir = pathlib.Path(res_dir)
            name = str(task.get('task_group')) + '_' + str(task.get('task_id'))
            data = task.get('res')
            data = get_msg_data(data, name=name)
            if data is None:
                data = b''
            write_bytes(res_dir / name, data)

    args = {'uids': set(uids), 'arg': arg}
    args.update(login)
    ## 注意队列maxsize太大比如80会将子进程卡主，目测应该是python的bug
    task_queue = mp.Queue(maxsize=50)  # 任务队列，投递任务
    result_queue = mp.Queue(maxsize=50)  # 结果队列，返回结果
    process_status_queue = mp.Queue(maxsize=50)  # 处理任务进程状态队列
    allocate_status_queue = mp.Queue(maxsize=50)  # 分配任务进程状态队列
    manager_queue = mp.Queue(maxsize=50)  # 处理任务进程管理，用于提示管理器结束子进程
    log_queue = mp.Queue()  # 任务记录队列，这里只记录了任务标识，没记录内容
    # 创建爬取任务分配进程
    print('assign email task')
    p1 = mp.Process(target=assign_mail_task,
                    kwargs={'task_queue': task_queue, 'log_queue': log_queue,
                            'status_queue': allocate_status_queue, 'args': args, 'in_parts': False})
    p1.start()
    # 创建邮件爬取任务处理进程
    print('create email process ==================')
    p0 = mp.Process(target=task_process_manager,
                    kwargs={'task_queue': task_queue, 'result_queue': result_queue,
                            'status_queue': process_status_queue, 'manager_queue': manager_queue, 'p_num': p_num,
                            'args': login})
    p0.start()
    # envents loop
    # wait log and task finish
    print('envents loop ===============')
    print('wait log and task finish, and get log and result')
    log_code = 1
    res_log = {}
    while True:
        if not log_queue.empty():
            log_code = 0
            # get log
            # log: {uid1: frag_total1,...}
            print('get log ===============')
            res_log = log_queue.get()
            if res_log == {}:
                raise ValueError('筛选任务为空，请检查邮箱和uid是否正确')
            if not os.path.exists(log_path):
                with open(log_path, 'w') as w:
                    for uid in res_log:
                        w.write(str(uid) + '\t' + str(res_log.get(uid)) + '\n')
        result_tmp = get_queue_all(result_queue)
        write_tmp(result_tmp=result_tmp, res_dir=res_dir)
        if log_code == 0 and get_status_code(process_status_queue) == 0:
            break
    # 关闭两个任务进程
    p1.join()
    p1.close()
    manager_queue.put(0)
    p0.join()
    p0.close()
    return res_log


def read_bytes(path):
    if os.path.exists(path):
        with open(path, 'rb') as f:
            data = f.read()
    else:
        data = None
    return data


def write_bytes(path, data: bytes):
    if data is None:
        data = b''
    with open(path, 'wb') as f:
        f.write(data)
    return 0


def get_status_code(status_queue: mp.Queue):
    status_tmp = get_queue_all(status_queue)
    if len(status_tmp):
        if status_tmp.pop() == 0:
            code = 0
        else:
            code = 1
    else:
        code = -1
    return code


def task_process_manager(task_queue: mp.Queue, result_queue: mp.Queue,
                         status_queue: mp.Queue, manager_queue: mp.Queue, args: dict, p_num: int = 4):
    process = []
    result_tmp = []
    timer = time.time()
    task_list = []
    # queue cache
    task_cache = 2
    result_cache = 2
    # create Email_task_process
    for v in range(p_num):
        process.append(Email_task_process(args=args, pid=v, task_cache=task_cache, result_cache=result_cache))
    # manager loop
    while True:
        if not manager_queue.empty():
            break
        if len(task_list) < 1000:
            task_list += get_queue_all(task_queue)
        status_tmp = len(task_list)
        for p in process:
            p.refresh_process_auto()
            p.refresh_timer_auto()
            if not p.is_task_full():
                put_queue_all(task_list, p.task_queue)
            result_tmp += get_queue_all(p.result_queue)
        put_queue_all(result_tmp, result_queue)
        if time.time() - timer > 1:
            while status_queue.full():
                get_queue_all(status_queue)
            status_queue.put(status_tmp)
            timer = time.time()
    # 结束子进程
    for p in process:
        p.close()
    return 0


class Email_task_process:
    def __init__(self, args: dict, pid: int, task_cache: int = 3, result_cache: int = 3) -> None:
        self.pid = pid
        self.args = args
        self.timeout = 30
        self.task_cache = task_cache
        self.result_cache = result_cache
        self.timer = time.time()
        self.task_queue = mp.Queue(maxsize=self.task_cache)
        self.result_queue = mp.Queue(maxsize=self.result_cache)
        self.process = mp.Process(target=email_task_deal_queue,
                                  kwargs={'task_queue': self.task_queue, 'result_queue': self.result_queue,
                                          'args': self.args, 'pid': self.pid})
        self.process.start()

    def refresh_timer_auto(self):
        if self.is_refresh_timer():
            self.refresh_timer()

    def is_refresh_timer(self):
        return self.is_result_full() or not self.is_task_full()

    def refresh_timer(self):
        self.timer = time.time()

    def is_timeout(self, timeout: Union[float, None] = None):
        if timeout is not None:
            self.timeout = timeout
        if time.time() - self.timer > self.timeout:
            return True
        else:
            return False

    def is_task_full(self):
        return self.task_queue.full()

    def is_result_full(self):
        return self.result_queue.full()

    def is_full(self):
        return self.is_result_full() or self.is_task_full()

    def is_refresh_process(self):
        return self.is_task_full() and not self.is_result_full() and self.is_timeout()

    def refresh_process_auto(self):
        if self.is_refresh_process():
            self.refresh_process()

    def refresh_process(self):
        print('refresh_process ' + str(self.pid))
        # close raw process
        self.process.kill()
        self.process.join(timeout=10)
        self.process.close()
        self.timer = time.time()
        # get queue data
        try:
            task_list = get_queue_all(self.task_queue)
        except:
            print('Email_task_process refresh_process  fail to get task_list ' + str(self.pid))
            task_list = []
        try:
            result_list = get_queue_all(self.result_queue)
        except:
            print('Email_task_process refresh_process  fail to get result_list ' + str(self.pid))
            result_list = []
        # close raw queue
        self.task_queue.close()
        self.result_queue.close()
        # create new queue
        self.task_queue = mp.Queue(maxsize=self.task_cache)
        self.result_queue = mp.Queue(maxsize=self.result_cache)
        put_queue_all(task_list, self.task_queue)
        put_queue_all(result_list, self.result_queue)
        # create new process
        self.process = mp.Process(target=email_task_deal_queue,
                                  kwargs={'task_queue': self.task_queue, 'result_queue': self.result_queue,
                                          'args': self.args, 'pid': self.pid})
        self.process.start()

    def close(self):
        # close raw queue
        self.task_queue.close()
        self.result_queue.close()
        # close raw process
        self.process.kill()
        self.process.join()
        self.process.close()


def put_queue_all(data: list, queue: mp.Queue):
    """向队列投递，投递完或者队列满了为止"""
    while True:
        if queue.full() or len(data) == 0:
            break
        else:
            task = data.pop()
            try:
                queue.put(task)
            except:
                data.append(task)
                break
    return 0


def get_queue_all(queue: mp.Queue):
    """获取队列全部内容，直到队列为空"""
    res = []
    while True:
        if queue.empty():
            break
        else:
            try:
                res.append(queue.get(timeout=3))
            except:
                break
    return res


def assign_mail_task(task_queue: mp.Queue, log_queue: mp.Queue, status_queue: mp.Queue,
                     args: dict, chunk_size: int = 10000, in_parts: bool = True):
    """
    邮件任务注册队列，向任务队列添加任务
    只有(RFC822)参数才能获取大小，分片下载邮件
    """
    uids = args.get('uids')
    base_arg = args.get('arg')
    if base_arg == '(RFC822)' and in_parts:
        assign_mail_task_in_parts(task_queue=task_queue, log_queue=log_queue, status_queue=status_queue, args=args,
                                  chunk_size=chunk_size)
    else:
        log_tmp = {}
        n = len(uids)
        timer = time.time()
        for uid in uids:
            # 推送状态
            n -= 1
            if time.time() - timer > 1:
                while status_queue.full():
                    get_queue_all(status_queue)
                status_queue.put(n)
                timer = time.time()
            # 生成任务
            task_group = uid
            args = {'uid': str(uid), 'arg': base_arg}
            task = {'task': 'fetch_mail_in_parts', 'task_id': 0, 'task_group': task_group, 'args': args}
            while task_queue.full():
                print('assign_mail_task is waiting')
                time.sleep(1)  # 可能会阻塞
            try:
                task_queue.put(task)  # 投递任务
            except:
                print('Fetch_mail_in_parts assign_task fail')
            log_tmp.update({uid: 0})
        # 推送状态
        while status_queue.full():
            get_queue_all(status_queue)
        status_queue.put(n)
        # 推送记录
        log_queue.put(log_tmp)
    # 结束进程
    status_queue.close()
    task_queue.close()
    log_queue.close()
    return 0


def assign_mail_task_in_parts(task_queue: mp.Queue, log_queue: mp.Queue, status_queue: mp.Queue,
                              args: dict, chunk_size: int = 5000):
    """
    邮件任务注册队列，向任务队列添加任务,分片下载邮件
    只有(RFC822)参数才能获取大小。
    """
    user = args.get('user')
    password = args.get('password')
    mailbox = args.get('mailbox')
    uids = args.get('uids')
    base_arg = args.get('arg')
    if base_arg == '(RFC822)':
        arg_size = re.sub('\)$', '.SIZE)', base_arg)
    else:
        raise ValueError('只有(RFC822)参数才能获取大小，分片下载邮件')
    email_server = Email(user=user, password=password, mailbox=mailbox)
    email_server.link_server()
    log_tmp = {}
    n = len(uids)
    timer = time.time()
    for uid in uids:
        # 推送状态
        n -= 1
        if time.time() - timer > 1:
            while status_queue.full():
                get_queue_all(status_queue)
            status_queue.put(n)
            timer = time.time()
        # 爬取邮件大小
        typ, msg_data = email_server.email_server.uid('fetch', str(uid), arg_size)  # 可能会阻塞
        if msg_data is None or msg_data[0] is None:
            print(str(uid) + ' None')
            continue
        msg_size = msg_data[0].split(b' ')[-1].decode()
        msg_size = int(re.sub('\)', '', msg_size))
        # 根据邮件大小拆分下载任务
        task_group = uid
        start_byte = 0
        end_byte = chunk_size - 1
        n = 0
        while True:
            if end_byte >= msg_size:
                end_byte = msg_size - 1
            # typ, body_data = email_server.uid('fetch', str(uid), f'(BODY.PEEK[]<{start_byte}.{end_byte}>)')
            args = {'uid': str(uid), 'arg': re.sub('\)$', '', base_arg) + f'<{start_byte}.{end_byte}>)'}
            # assign_task
            task = {'task': 'fetch_mail_in_parts', 'task_id': n, 'task_group': task_group, 'args': args}
            while task_queue.full():
                time.sleep(1)  # 可能会阻塞
            try:
                task_queue.put(task)  # 投递任务
            except:
                print('Fetch_mail_in_parts assign_task fail')
                continue
            if end_byte == msg_size - 1:
                break
            start_byte = end_byte + 1
            end_byte = start_byte + chunk_size - 1
            n += 1
        log_tmp.update({uid: n - 1})
    # 推送记录
    log_queue.put(log_tmp)
    # 推送状态
    while status_queue.full():
        get_queue_all(status_queue)
    status_queue.put(n)
    # 结束进程
    status_queue.close()
    task_queue.close()
    log_queue.close()
    return 0


def email_task_deal_queue(task_queue: mp.Queue, result_queue: mp.Queue, args: dict, pid=None):
    """邮件任务处理队列，在任务队列中获取任务并执行"""
    user = args.get('user')
    password = args.get('password')
    mailbox = args.get('mailbox')
    email_server = Email(user=user, password=password, mailbox=mailbox)
    email_server.link_server()
    while True:
        if task_queue.empty() or result_queue.full():
            print("process " + str(pid) + " is waiting")
            time.sleep(1)
            continue
        else:
            try:
                task = task_queue.get(timeout=3)
            except:
                print('email_task get task fail ' + str(pid))
                continue
            task_id = task.get('task_id')
            task_group = task.get('task_group')
            print("process " + str(pid) + " deal, id: " + str(task_id) + ' group: ' + str(task_group))
            task_args = task.get('args')
            email_server.check_server()
            res = email_server.fetch_by_uid(uid=task_args.get('uid'), arg=task_args.get('arg'))
            task.update({'res': res})
            result_queue.put(task)
            print("process " + str(pid) + " deal, id: " + str(task_id) + ' group: ' + str(task_group) + 'done')


class Email:
    def __init__(self,
                 user: Union[str, None] = None,
                 password: Union[str, None] = None,
                 mailbox: Union[str, None] = None
                 ) -> None:
        self.user = user
        self.password = password
        self.mailbox = mailbox

    def link_server(self):
        """建立与邮箱服务器的连接"""
        user = self.user
        password = self.password
        mailbox = self.mailbox
        if user is None or password is None or mailbox is None:
            raise ValueError('login : user is None or password is None')
        email_server = IMAP4_SSL(host='imap.exmail.qq.com', port=993)
        email_server.login(user=user, password=password)
        email_server.select(mailbox, readonly=True)
        # email_server.select(readonly=True)
        self.email_server = email_server

    def check_server(self):
        """检查连接"""
        try:
            self.email_server.noop()
        except:
            self.link_server()

    def fetch_by_uid(self, uid, arg):
        typ, msg_data = self.email_server.uid('fetch', str(uid), arg)
        return msg_data

    def uid(self, **key):
        return self.email_server.uid(**key)

    def select(self, **key):
        self.email_server.select(**key)


def merge_result(result: dict, log: dict):
    """
    合并分片下载邮件结果
    result, 嵌套字典:{uid1: {frag1: body1,...},...}
    log: {uid1: frag_total1,...}
    return {uid: body}
    """
    res = {}
    for uid in log:
        frag_total = result.get(uid)
        body = b''
        if frag_total is not None:
            for frag in range(len(frag_total)):
                body_frag = frag_total.get(frag)
                if isinstance(body_frag, bytes):
                    body += body_frag
                elif isinstance(body_frag, list):
                    if isinstance(body_frag[0], tuple):
                        body += body_frag[0][1]
                else:
                    print(str(uid) + '-' + str(frag) + ' missing')
        else:
            print(str(uid) + ' missing')
        res.update({uid: body})
    return res


def decode_header(msg, con: Union[str, None] = None, join: Union[str, None] = None, pattern: Union[str, None] = None):
    """对email.message.Message对象的元素进行字符集转化，空元素和转化失败的都为None"""
    # email.header.decode_header 尝试识别字符集类型，可能会出现识别错误
    if con is None and msg is not None:
        res = [try_parse(x) for x in email.header.decode_header(msg)]
    elif msg is None or msg[con] is None:
        res = [None]
    else:
        res = [try_parse(x) for x in email.header.decode_header(msg[con])]
    res = [x for x in res if x is not None]
    if pattern is not None:
        res = [x for x in res if re.search(pattern, x)]
    if join is not None:
        res = join.join([x for x in res if x is not None])
    return res


def try_parse(data: tuple, err: bool = False):
    """对email.header.decode_header识别的字符尝试按照识别字符集类型转化，失败则通用转化"""
    if len(data) == 2 and data[1] is not None:
        res = is_decode(data=data[0], code=data[1])
        if res is None:
            res = try_decode(data[0], err=err)
    else:
        res = try_decode(data[0], err=err)
    return res


def is_decode(data, code: Union[str, None] = None, err: bool = False):
    """尝试按照指定编码格式进行字符集转化"""
    try:
        if code is None:
            res = data.decode()
        else:
            res = data.decode(code)
    except:
        if err:
            raise UnicodeDecodeError('decode fail')
        else:
            res = None
    return res


def try_decode(data, err=False):
    """通用字符集转化，依次尝试多种编码"""
    if isinstance(data, str):
        res = data
    else:
        res = is_decode(data, code=None)
        if res is None:
            res = is_decode(data, code='utf-8')
        if res is None:
            res = is_decode(data, code='gb2312')
        if res is None:
            res = is_decode(data, code='gbk')
        if res is None and err:
            raise UnicodeDecodeError('decode fail')
    return res


def add_timezone(dt: datetime.datetime):
    """没有时区信息，添加时区信息"""
    if dt.tzinfo is None or dt.tzinfo.utcoffset(dt) is None:
        tz = pytz.timezone('Asia/Shanghai')
        dt = tz.localize(dt)
    return dt


def remove_timezone(dt: Union[datetime.datetime, time.struct_time]):
    """
    去掉datetime对象的时区信息
    """
    if isinstance(dt, time.struct_time):
        dt = datetime.datetime.fromtimestamp(time.mktime(dt))
    if dt.tzinfo is not None:
        dt = dt.replace(tzinfo=None)
    return dt


def message_parsing(message):
    """parsing messages into a readable string"""
    messages = email.message_from_bytes(message)
    subject = decode_header(messages, "Subject", join=',', pattern=None)
    print(subject)
    soup_list = []
    file_list = []
    for data in messages.walk():
        if data.get_content_type() == 'text/plain':
            try:
                soup = BeautifulSoup(base64.b64decode(try_decode(data.get_payload(decode=True), err=True)),
                                     'html.parser').text
            except:
                try:
                    soup = BeautifulSoup(base64.b64decode(data.get_payload()), 'html.parser').text
                except:
                    soup = data.get_payload()
            soup = soup.replace('\r\n', '\n').replace('\xa0', ' ')
            soup_list.append(soup)
            a = re.search("从腾讯企业邮箱发来的超大附件\r\n\r\n(.+?)\s*\(.+", soup)
            if a is not None:
                file_list.append(a.group(1))
        else:
            if re.search('application', data.get_content_type()) and data.get('Content-Disposition') is not None:
                filename = decode_header(data.get_filename())
                if filename is not None:
                    name = ",".join(filename)
                    file_list.append(name)
    soup_list = [x for x in soup_list if x is not None]
    file_list = [x for x in file_list if x is not None]
    return soup_list, file_list


if __name__ == '__main__':
    main_task()
    # con = [{'name': '数据', 'SINCE': '30-Nov-2021', 'SUBJECT': '数据'}, 
    #     {'name': '报告', 'SINCE': '30-Nov-2022', 'SUBJECT': '报告'}]
    # con1 = [{'name': '数据', 'SINCE': '30-Nov-2022', 'SUBJECT': '数据'}]
    # con2 = [{'name': '报告', 'SINCE': '30-Nov-2022', 'SUBJECT': '报告'}]
    # print('main_task1')
    # main_task(user='xing.fang@oebiotech.com', password='F2384x', mailbox='"Sent Messages"', con=con, load=True)
    # print('main_task2')
    # main_task(user='xing.fang@oebiotech.com', password='F2384x', mailbox='INBOX',con=con, load=True)
    # print('main_task3')
    # main_task(user='hui.zhu@oebiotech.com', password='6hLTVRchnDMhDzW', mailbox='INBOX', con=con2, load=True)
    # print('main_task4')
    # main_task(user='hui.zhu@oebiotech.com', password='6hLTVRchnDMhDzW', mailbox='"Sent Messages"', con=con1, load=True)
    # print('main_task5')
    # main_task(user='xinyue.zhang@oebiotech.com', password='OEzhangxinyue99', mailbox='"Sent Messages"', con=con2, load=True)
    # print('main_task6')
    # main_task(user='fang.qiao@oebiotech.com', password='6450676Qiao', mailbox='"Sent Messages"', con=con2, load=True)
    # print('main_task7')
    # main_task(user='qinqin.zhu@oebiotech.com', password='zm941213QQ', mailbox='"Sent Messages"', con=con2, load=True)
    # print('main_task8')
    # main_task(user='yalan.zhao@oebiotech.com', password='Quan19880818', mailbox='"Sent Messages"', con=con2, load=True)
    # print('main_task9')
    # main_task(user='lishu.chen@oebiotech.com', password='1999ZYcl', mailbox='"Sent Messages"', con=con1, load=True)
    
