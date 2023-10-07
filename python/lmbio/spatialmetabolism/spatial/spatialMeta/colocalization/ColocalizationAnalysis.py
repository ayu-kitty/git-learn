# File: ColocalizationAnalysis.py
# Author: jiang tao
# Time: 2023/3/7 11:20
import logging
import re
from typing import List, Union

import attr
from pathlib import Path

from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
import pandas as pd

from .picMatch import rotateHEMatchMsi
from .resnet import exportResetFeatures
from .rutils import (base, lmbio, cardinal, ebi_image,
                     to_rdf, fig_deal, get_fig_data_for_resnet,
                     spatial_rawdata_metabolite_correlation, spatial_resnetdata_figure_correlation,
                     NULL)

logger = logging.getLogger("SpatialMeta")
MODULE = "空代共定位分析: "


@attr.s(auto_attribs=True)
class ColocalizationAnalysis:
    output_folder: str
    imzml_path: Union[str, List[str]]
    he_path: Union[str, List[str]]

    @staticmethod
    def parse_sample_info(imzml_path: str, he_path: str):
        # eg. L-1-neg.imzML
        imzml = Path(imzml_path).name
        # eg. L-1
        sample_name = re.sub(r"-\w{3}\.imzML", "", imzml)
        # Ion mode
        match = re.search(r"-neg\.", imzml, flags=re.I)
        ion_mode = "neg" if match else "pos"
        # eg. LPDL-1
        fig_name = Path(he_path).stem
        return sample_name, fig_name, ion_mode

    @staticmethod
    def rdf_pandas(rdf):
        with (ro.default_converter + pandas2ri.converter).context():
            df = ro.conversion.get_conversion().rpy2py(rdf)
        return df

    @staticmethod
    def pandas_rdf(pd_df):
        with (ro.default_converter + pandas2ri.converter).context():
            rdf = ro.conversion.get_conversion().py2rpy(pd_df)
        return rdf

    def getCoordinates(self, imzml_path: str) -> pd.DataFrame:
        """从imzml中读取成像图的坐标数据"""
        logger.info(MODULE + "从imzml中读取成像图的坐标数据...")
        mse = cardinal.readMSIData(imzml_path)
        # RS4
        coord_data = cardinal.coord(mse)
        # rdf
        coord_data = to_rdf(coord_data)
        return self.rdf_pandas(coord_data)

    def heMatch(self, imzml_path: str, he_path: str, data_msi: pd.DataFrame,
                image_type: str = 'png'):
        """
        @Step1: 旋转`HE染色图`与`代谢物成像图`[切片]匹配
        :param imzml_path: imzml abs path
        :param data_msi: 代谢物成像图 df
        :param he_path: HE染色图 abs path
        :param ion_mode: ion mode
        :param image_type: image type
        :return 自动匹配图像路径 fig_file
                在out_path生成最佳匹配的HE，或者排名前十
        """
        logger.info(MODULE + "自动旋转`HE染色图`与`代谢物成像图`匹配...")
        out_path = self.output_folder
        sample_name, fig_name, ion_mode = self.parse_sample_info(imzml_path, he_path)
        out_f = Path(out_path) / "sample" / "fig" / "match" / sample_name / ion_mode / fig_name
        if not out_f.exists():
            out_f.mkdir(parents=True, exist_ok=True)
        rotateHEMatchMsi(data_msi, he_path, ion_mode, out_f)

        # imzmlimage 额外导出一张图
        fig_file = str(out_f / f"cropped_BestMatch_{ion_mode}.png")
        lmbio.imzmlimage(filename=imzml_path,
                         savepath=str(Path(out_path) / "sample" / "fig"),
                         figfile=fig_file,
                         mapname=sample_name + "-" + ion_mode + "-" + fig_name,
                         imagetype=image_type)
        print("自动匹配图像路径:", fig_file)
        return fig_file

    def heMatchManual(self, imzml_path: str, he_path: str, fig_file: str,
                      image_type: str = 'png',
                      ps=5,
                      tv=0.1,
                      bg_alpha=0,
                      fig_alpha=0.2,
                      pixel_thresh=0.85, num=5,
                      ymax_adj=0, ymin_adj=0,
                      xmax_adj=0, xmin_adj=0,
                      ymax_adj_=0, ymin_adj_=0,
                      xmax_adj_=0, xmin_adj_=0,
                      interaction=False,
                      rotate_angle=0,
                      flip=False,
                      flop=False,
                      resize=False,
                      fig_width=1000,
                      fig_height=1000
                      ):
        """手动调整HE染色图"""
        logger.info(MODULE + "手动调整`HE染色图`与`代谢物成像图`匹配...")
        out_path = self.output_folder
        sample_name, fig_name, ion_mode = self.parse_sample_info(imzml_path, he_path)
        out_dir = Path(out_path) / "sample" / "fig" / "process" / sample_name / ion_mode / fig_name
        if not out_dir.exists():
            out_dir.mkdir(parents=True, exist_ok=True)
        # fig_data rpy2.robjects.vectors.ListVector
        fig_data = fig_deal(file=str(Path(fig_file)),
                            savepath=str(out_dir),
                            ps=ps,
                            tv=tv,
                            bg_alpha=bg_alpha,
                            fig_alpha=fig_alpha,
                            pixel_thresh=pixel_thresh, num=num,
                            ymax_adj=ymax_adj, ymin_adj=ymin_adj,
                            xmax_adj=xmax_adj, xmin_adj=xmin_adj,
                            ymax_adj_=ymax_adj_, ymin_adj_=ymin_adj_,
                            xmax_adj_=xmax_adj_, xmin_adj_=xmin_adj_,
                            interaction=interaction,
                            rotate_angle=rotate_angle,
                            flip=flip,
                            flop=flop,
                            resize=resize,
                            width=fig_width,
                            height=fig_height)

        fig_cut_file = str(
            Path(out_path) / "sample" / "fig" / (sample_name + "-" + ion_mode + "-" + fig_name + "-fig.png"))
        fig_data.rx2['figfile'] = fig_cut_file
        ebi_image.writeImage(x=fig_data.rx2['figcut'], files=fig_cut_file, type=image_type)

        rds_file = str(Path(out_path) / "sample" / "fig" / (sample_name + "-" + ion_mode + "-" + fig_name + ".rds"))
        fig_data.rx2['rdsfile'] = rds_file
        fig_data.rx2['imzmlfile'] = imzml_path
        base.saveRDS(fig_data, file=rds_file)

        lmbio.imzmlimage(filename=imzml_path,
                         savepath=str(Path(out_path) / "sample" / "fig"),
                         figfile=fig_cut_file,
                         mapname=sample_name + "-" + ion_mode + "-" + fig_name + "-manual",
                         imagetype=image_type)

    def get_fig_data(self, imzml_path: str, he_path: str):
        logger.info(MODULE + "提取fig data for resnet...")
        out_path = self.output_folder
        sample_name, fig_name, ion_mode = self.parse_sample_info(imzml_path, he_path)
        rds_file = str(Path(out_path) / "sample" / "fig" / (sample_name + "-" + ion_mode + "-" + fig_name + ".rds"))
        # R脚本 get_fig_data_for_resnet
        # return：R list ; this step error in windows
        fig_data = get_fig_data_for_resnet(imzml_path, mode=ion_mode, samplename=sample_name,
                                           figrds=ro.StrVector((rds_file,)),
                                           figname=ro.StrVector((fig_name,))
                                           )
        return fig_data

    def get_resnet_data(self, imzml_path: str, he_path: str, fig_data: ro.ListVector, parallel: bool = True) -> pd.DataFrame:
        logger.info(MODULE + "从图片提取特征 resnet-2048...")
        out_path = self.output_folder
        sample_name, fig_name, ion_mode = self.parse_sample_info(imzml_path, he_path)
        # out_path
        out_resnet = Path(out_path) / "sample" / "fig" / "resnet" / sample_name / ion_mode / fig_name
        if not out_resnet.exists():
            out_resnet.mkdir(parents=True, exist_ok=True)
        # r,g,b 随机挑了一个 r
        data = self.rdf_pandas(fig_data.rx2['anadata'].rx2['r'])
        # TODO: 直接返回数据，不写出 resnet_data
        resnet_data = exportResetFeatures(data, str(out_resnet), draw=True, parallel=parallel)
        return resnet_data


    def correlation_raw_metabolite(self, imzml_path: str, he_path: str, fig_data: ro.ListVector,
                                   selectmz=NULL,
                                   corfilter=0.7,
                                   corfiltertype="+",
                                   pfilter=0.05,
                                   adjust="BH",
                                   cormethod="pearson",
                                   use="pairwise",
                                   xname="mz-x",
                                   yname="mz-y"):
        logger.info(MODULE + "绘制代谢物相关性FROM RAW DATA...")
        out_path = self.output_folder
        sample_name, fig_name, ion_mode = self.parse_sample_info(imzml_path, he_path)

        out_meta = Path(out_path) / "sample" / "fig" / "result" / "meta"
        if not out_meta.exists():
            out_meta.mkdir(parents=True, exist_ok=True)

        if isinstance(selectmz, list):
            selectmz = ro.StrVector(selectmz)
        spatial_rawdata_metabolite_correlation(
            imzml_path,
            df=fig_data.rx2['spectradata'],
            mode=ion_mode,
            samplename=sample_name,
            selectmz=selectmz,
            savepath=str(out_meta),
            corfilter=corfilter,
            corfiltertype=corfiltertype,
            pfilter=pfilter,
            adjust=adjust,
            cormethod=cormethod,
            use=use,
            xname=xname,
            yname=yname)

    def correlation_resnet(self, imzml_path: str, he_path: str, resnet_data: pd.DataFrame,
                                   selectmz=NULL,
                                   corfilter=0.7,
                                   corfiltertype="+",
                                   pfilter=0.05,
                                   adjust="BH",
                                   cormethod="pearson",
                                   use="pairwise",
                                   xname="mz-x",
                                   yname="mz-y"):

        logger.info(MODULE + "绘制HE与代谢物相关性FROM RESNET DATA...")
        out_path = self.output_folder
        sample_name, fig_name, ion_mode = self.parse_sample_info(imzml_path, he_path)

        out_meta = Path(out_path) / "sample" / "fig" / "result" / "resnet"
        if not out_meta.exists():
            out_meta.mkdir(parents=True, exist_ok=True)
        if isinstance(selectmz, list):
            selectmz = ro.StrVector(selectmz)

        resnet_data = self.pandas_rdf(resnet_data)
        spatial_resnetdata_figure_correlation(
            imzml_path,
            resnet_data,
            mode=ion_mode,
            samplename=sample_name,
            selectmz=selectmz,
            savepath=str(out_meta),
            corfilter=corfilter,
            corfiltertype=corfiltertype,
            pfilter=pfilter,
            adjust=adjust,
            cormethod=cormethod,
            use=use,
            xname=xname,
            yname=yname)

    def run_pipeline(self, selectmz=NULL,
                           corfilter=0.7,
                           corfiltertype="+",
                           pfilter=0.05,
                           adjust="BH",
                           cormethod="pearson",
                           use="pairwise",
                           xname="mz-x",
                           yname="mz-y"):
        """空代共定位分析的完整 pipline"""
        logger.info(MODULE + "run_pipeline...")
        imzml_path, he_path = self.imzml_path, self.he_path
        if not isinstance(self.imzml_path, list):
            imzml_path, he_path = [self.imzml_path], [self.he_path]
        for imzml, he in zip(imzml_path, he_path):
            logger.info(MODULE + f"正在处理{imzml}...")
            df = self.getCoordinates(imzml)
            fig_file = self.heMatch(imzml, he, df)
            self.heMatchManual(imzml, he, fig_file)
            fig_data = self.get_fig_data(imzml, he)
            self.correlation_raw_metabolite(imzml, he, fig_data,
                                            selectmz=selectmz,
                                            corfilter=corfilter,
                                            corfiltertype=corfiltertype,
                                            pfilter=pfilter,
                                            adjust=adjust,
                                            cormethod=cormethod,
                                            use=use,
                                            xname=xname,
                                            yname=yname)
            resnet_data = self.get_resnet_data(imzml, he, fig_data)
            self.correlation_resnet(imzml, he, resnet_data,
                                    selectmz=selectmz,
                                    corfilter=corfilter,
                                    corfiltertype=corfiltertype,
                                    pfilter=pfilter,
                                    adjust=adjust,
                                    cormethod=cormethod,
                                    use=use,
                                    xname=xname,
                                    yname=yname)
        logger.info(MODULE + "run_pipeline 分析完成!")

    def manual_pipeline(self,  selectmz=NULL,
                               corfilter=0.7,
                               corfiltertype="+",
                               pfilter=0.05,
                               adjust="BH",
                               cormethod="pearson",
                               use="pairwise",
                               xname="mz-x",
                               yname="mz-y"):
        """手动调整HE染色图后， 直接运行这一步完成后续分析"""
        logger.info(MODULE + " 手动调整HE染色图后 - run manual pipeline...")
        imzml_path, he_path = self.imzml_path, self.he_path
        if not isinstance(self.imzml_path, list):
            imzml_path, he_path = [self.imzml_path], [self.he_path]
        for imzml, he in zip(imzml_path, he_path):
            logger.info(MODULE + f"正在处理{imzml}...")
            fig_data = self.get_fig_data(imzml, he)
            self.correlation_raw_metabolite(imzml, he, fig_data,
                                            selectmz=selectmz,
                                            corfilter=corfilter,
                                            corfiltertype=corfiltertype,
                                            pfilter=pfilter,
                                            adjust=adjust,
                                            cormethod=cormethod,
                                            use=use,
                                            xname=xname,
                                            yname=yname)
            resnet_data = self.get_resnet_data(imzml, he, fig_data)
            self.correlation_resnet(imzml, he, resnet_data,
                                    selectmz=selectmz,
                                    corfilter=corfilter,
                                    corfiltertype=corfiltertype,
                                    pfilter=pfilter,
                                    adjust=adjust,
                                    cormethod=cormethod,
                                    use=use,
                                    xname=xname,
                                    yname=yname)
        logger.info(MODULE + "run manual pipeline 分析完成!")


if __name__ == '__main__':
    pass
    # ca = ColocalizationAnalysis(output_folder="/home/jiangt/pyCharm/spatial/tests/output")
    # imzml_path = "/home/jiangt/pyCharm/spatial/tests/imzml/L-1-neg.imzML"
    # he_path = "/home/jiangt/pyCharm/spatial/tests/imzml/LPDL-1.jpg"
    # df = ca.getCoordinates(imzml_path)
    # fig_file = ca.heMatch(imzml_path, he_path, df)
    # ca.heMatchManual(imzml_path, he_path, fig_file)
    # fig_data = ca.get_fig_data(imzml_path, he_path)
    # ca.correlation_raw_metabolite(imzml_path, he_path, fig_data)
    # resnet_data = ca.get_resnet_data(imzml_path, he_path, fig_data)
    # ca.correlation_resnet(imzml_path, he_path, resnet_data)
