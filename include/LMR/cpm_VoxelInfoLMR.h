/*
###################################################################################
#
# CPMlib - Computational space Partitioning Management library
#
# Copyright (c) 2012-2014 Institute of Industrial Science (IIS), The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2014-2016 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
 */

/**
 * @file   cpm_VoxelInfoLMR.h
 * LMR用のVOXEL空間情報クラスのヘッダーファイル
 * @date   2015/03/27
 */

#ifndef _CPM_VOXELINFO_LMR_H_
#define _CPM_VOXELINFO_LMR_H_

#include "cpm_VoxelInfo.h"
//#include "BCM/BCMOctree.h"
//#include "BCM/BCMFileCommon.h"
#include "BCMOctree.h"
#include "BCMFileCommon.h"
#include "cpm_TextParserDomainLMR.h"
#include "cpm_DefineLMR.h"

///** リーフ毎のVOXEL空間情報管理マップ */
class cpm_VoxelInfoLMR;
typedef std::map<int, cpm_VoxelInfoLMR*> LeafMap; //map<leafID,VoxelInfo*>

/** LMR用のVOXEL空間情報管理クラス
 */
class cpm_VoxelInfoLMR : public cpm_VoxelInfo
{
friend class cpm_BaseParaManager;
friend class cpm_ParaManagerLMR;
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  void debugPrint()
  {
    const int *gv = m_globalDomainInfo.GetVoxNum();
    const int *dv = GetDivNum();
    const int *ps = GetDivPos();
    const int *lv = m_localDomainInfo.GetVoxNum();
    const int *hd = m_voxelHeadIndex;
    const int *tl = m_voxelTailIndex;
    std::cout << "*** leaf  = " << m_leafID << std::endl;
    std::cout << "    rank  = " << m_rankNo << std::endl;
    std::cout << "    level = " << m_node->getLevel() << std::endl;
    std::cout << "    gvox  = " << gv[0] << "," << gv[1] << "," << gv[2] << std::endl;
    std::cout << "    div   = " << dv[0] << "," << dv[1] << "," << dv[2] << std::endl;
    std::cout << "    pos   = " << ps[0] << "," << ps[1] << "," << ps[2] << std::endl;
    std::cout << "    lvox  = " << lv[0] << "," << lv[1] << "," << lv[2] << std::endl;
    std::cout << "    head  = " << hd[0] << "," << hd[1] << "," << hd[2] << std::endl;
    std::cout << "    tail  = " << tl[0] << "," << tl[1] << "," << tl[2] << std::endl;
    std::cout << "    -X Rn = " << m_neighborRankID_LMR[X_MINUS][0] << "," << m_neighborRankID_LMR[X_MINUS][1] << "," << m_neighborRankID_LMR[X_MINUS][2] << "," << m_neighborRankID_LMR[X_MINUS][3] << std::endl;
    std::cout << "    +X Rn = " << m_neighborRankID_LMR[X_PLUS ][0] << "," << m_neighborRankID_LMR[X_PLUS ][1] << "," << m_neighborRankID_LMR[X_PLUS ][2] << "," << m_neighborRankID_LMR[X_PLUS ][3] << std::endl;
    std::cout << "    -Y Rn = " << m_neighborRankID_LMR[Y_MINUS][0] << "," << m_neighborRankID_LMR[Y_MINUS][1] << "," << m_neighborRankID_LMR[Y_MINUS][2] << "," << m_neighborRankID_LMR[Y_MINUS][3] << std::endl;
    std::cout << "    +Y Rn = " << m_neighborRankID_LMR[Y_PLUS ][0] << "," << m_neighborRankID_LMR[Y_PLUS ][1] << "," << m_neighborRankID_LMR[Y_PLUS ][2] << "," << m_neighborRankID_LMR[Y_PLUS ][3] << std::endl;
    std::cout << "    -Z Rn = " << m_neighborRankID_LMR[Z_MINUS][0] << "," << m_neighborRankID_LMR[Z_MINUS][1] << "," << m_neighborRankID_LMR[Z_MINUS][2] << "," << m_neighborRankID_LMR[Z_MINUS][3] << std::endl;
    std::cout << "    +Z Rn = " << m_neighborRankID_LMR[Z_PLUS ][0] << "," << m_neighborRankID_LMR[Z_PLUS ][1] << "," << m_neighborRankID_LMR[Z_PLUS ][2] << "," << m_neighborRankID_LMR[Z_PLUS ][3] << std::endl;
    std::cout << "    -X Rp = " << m_periodicRankID_LMR[X_MINUS][0] << "," << m_periodicRankID_LMR[X_MINUS][1] << "," << m_periodicRankID_LMR[X_MINUS][2] << "," << m_periodicRankID_LMR[X_MINUS][3] << std::endl;
    std::cout << "    +X Rp = " << m_periodicRankID_LMR[X_PLUS ][0] << "," << m_periodicRankID_LMR[X_PLUS ][1] << "," << m_periodicRankID_LMR[X_PLUS ][2] << "," << m_periodicRankID_LMR[X_PLUS ][3] << std::endl;
    std::cout << "    -Y Rp = " << m_periodicRankID_LMR[Y_MINUS][0] << "," << m_periodicRankID_LMR[Y_MINUS][1] << "," << m_periodicRankID_LMR[Y_MINUS][2] << "," << m_periodicRankID_LMR[Y_MINUS][3] << std::endl;
    std::cout << "    +Y Rp = " << m_periodicRankID_LMR[Y_PLUS ][0] << "," << m_periodicRankID_LMR[Y_PLUS ][1] << "," << m_periodicRankID_LMR[Y_PLUS ][2] << "," << m_periodicRankID_LMR[Y_PLUS ][3] << std::endl;
    std::cout << "    -Z Rp = " << m_periodicRankID_LMR[Z_MINUS][0] << "," << m_periodicRankID_LMR[Z_MINUS][1] << "," << m_periodicRankID_LMR[Z_MINUS][2] << "," << m_periodicRankID_LMR[Z_MINUS][3] << std::endl;
    std::cout << "    +Z Rp = " << m_periodicRankID_LMR[Z_PLUS ][0] << "," << m_periodicRankID_LMR[Z_PLUS ][1] << "," << m_periodicRankID_LMR[Z_PLUS ][2] << "," << m_periodicRankID_LMR[Z_PLUS ][3] << std::endl;
    std::cout << "    -X Ln = " << m_neighborLeafID_LMR[X_MINUS][0] << "," << m_neighborLeafID_LMR[X_MINUS][1] << "," << m_neighborLeafID_LMR[X_MINUS][2] << "," << m_neighborLeafID_LMR[X_MINUS][3] << std::endl;
    std::cout << "    +X Ln = " << m_neighborLeafID_LMR[X_PLUS ][0] << "," << m_neighborLeafID_LMR[X_PLUS ][1] << "," << m_neighborLeafID_LMR[X_PLUS ][2] << "," << m_neighborLeafID_LMR[X_PLUS ][3] << std::endl;
    std::cout << "    -Y Ln = " << m_neighborLeafID_LMR[Y_MINUS][0] << "," << m_neighborLeafID_LMR[Y_MINUS][1] << "," << m_neighborLeafID_LMR[Y_MINUS][2] << "," << m_neighborLeafID_LMR[Y_MINUS][3] << std::endl;
    std::cout << "    +Y Ln = " << m_neighborLeafID_LMR[Y_PLUS ][0] << "," << m_neighborLeafID_LMR[Y_PLUS ][1] << "," << m_neighborLeafID_LMR[Y_PLUS ][2] << "," << m_neighborLeafID_LMR[Y_PLUS ][3] << std::endl;
    std::cout << "    -Z Ln = " << m_neighborLeafID_LMR[Z_MINUS][0] << "," << m_neighborLeafID_LMR[Z_MINUS][1] << "," << m_neighborLeafID_LMR[Z_MINUS][2] << "," << m_neighborLeafID_LMR[Z_MINUS][3] << std::endl;
    std::cout << "    +Z Ln = " << m_neighborLeafID_LMR[Z_PLUS ][0] << "," << m_neighborLeafID_LMR[Z_PLUS ][1] << "," << m_neighborLeafID_LMR[Z_PLUS ][2] << "," << m_neighborLeafID_LMR[Z_PLUS ][3] << std::endl;
    std::cout << "    -X Lp = " << m_periodicLeafID_LMR[X_MINUS][0] << "," << m_periodicLeafID_LMR[X_MINUS][1] << "," << m_periodicLeafID_LMR[X_MINUS][2] << "," << m_periodicLeafID_LMR[X_MINUS][3] << std::endl;
    std::cout << "    +X Lp = " << m_periodicLeafID_LMR[X_PLUS ][0] << "," << m_periodicLeafID_LMR[X_PLUS ][1] << "," << m_periodicLeafID_LMR[X_PLUS ][2] << "," << m_periodicLeafID_LMR[X_PLUS ][3] << std::endl;
    std::cout << "    -Y Lp = " << m_periodicLeafID_LMR[Y_MINUS][0] << "," << m_periodicLeafID_LMR[Y_MINUS][1] << "," << m_periodicLeafID_LMR[Y_MINUS][2] << "," << m_periodicLeafID_LMR[Y_MINUS][3] << std::endl;
    std::cout << "    +Y Lp = " << m_periodicLeafID_LMR[Y_PLUS ][0] << "," << m_periodicLeafID_LMR[Y_PLUS ][1] << "," << m_periodicLeafID_LMR[Y_PLUS ][2] << "," << m_periodicLeafID_LMR[Y_PLUS ][3] << std::endl;
    std::cout << "    -Z Lp = " << m_periodicLeafID_LMR[Z_MINUS][0] << "," << m_periodicLeafID_LMR[Z_MINUS][1] << "," << m_periodicLeafID_LMR[Z_MINUS][2] << "," << m_periodicLeafID_LMR[Z_MINUS][3] << std::endl;
    std::cout << "    +Z Lp = " << m_periodicLeafID_LMR[Z_PLUS ][0] << "," << m_periodicLeafID_LMR[Z_PLUS ][1] << "," << m_periodicLeafID_LMR[Z_PLUS ][2] << "," << m_periodicLeafID_LMR[Z_PLUS ][3] << std::endl;
#if 0
    std::cout << "" <<  << std::endl;
    std::cout << "" <<  << std::endl;
    std::cout << "" <<  << std::endl;
    std::cout << "" <<  << std::endl;
    std::cout << "" <<  << std::endl;
#endif
  }


protected:
  /** コンストラクタ */
  cpm_VoxelInfoLMR();

  /** デストラクタ */
  virtual ~cpm_VoxelInfoLMR();


////// 領域分割情報生成関係 //////

  /** CPM領域分割情報の生成
   *  - MPI_COMM_WORLDを使用した領域を生成する。
   *
   *  @param[in]  comm    MPIコミュニケータ
   *  @param[in]  treeFile  領域情報ファイル
   *  @param[out] leafMap リーフごとのVoxel空間情報マップ
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  static cpm_ErrorCode Init( MPI_Comm comm, std::string treeFile, LeafMap &leafMap );

  /** 木情報ファイルの読み込み
   *  @param[in]  octFile   木情報ファイル
   *  @param[out] header    ヘッダー情報
   *  @param[out] pedigrees ぺディグリー情報
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  static
  cpm_ErrorCode LoadOctreeFile( std::string octFile, BCMFileIO::OctHeader &header, std::vector<Pedigree> &pedigrees );

  /** 木情報ファイルのヘッダー読み込み
   *  - ヘッダーのみを読み込み、ファイルをクローズする
   *
   *  @param[in]  octFile 木情報ファイル
   *  @param[out] header  ヘッダー情報
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  static
  cpm_ErrorCode LoadOctreeHeader( std::string octFile, BCMFileIO::OctHeader &header );

  /** 木情報ファイルのヘッダー読み込み
   *  @param[in]  fp         木情報ファイルポインタ
   *  @param[out] header     ヘッダー情報
   *  @param[out] isNeedSwap エンディアン変換フラグ(true:要変換)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  static
  cpm_ErrorCode LoadOctreeHeader( FILE *fp, BCMFileIO::OctHeader &header, bool &isNeedSwap );

  /** 各ランクの担当リーフマップを取得
   *  @param[in] nRank   並列数
   *  @param[in] numLeaf 全リーフ数
   *  @return リーフマップ(map<leafID,rankNo>)
   */
  static
  std::map<int,int> GetLeafIDMap(int nRank, int numLeaf);

  /** グローバルの領域情報をセット
   *  @param[in] dInfo 領域情報
   */
  void SetGlobaliDomainInfo( S_OCT_DOMAIN_INFO &dInfo );

  /** ローカルの領域情報をセット
   *  @param[in] dInfo 領域情報
   */
  void SetLocalDomainInfo( S_OCT_DOMAIN_INFO &dInfo );

  /** 隣接情報の取得
   */
  void SetNeighborInfo(const std::map<int,int> &leafIDmap);

  /** 木情報ファイルからリーフ数を取得する
   *  @param[in] treefile  木情報ファイル
   *  @return    リーフ数
   */
  static
  int GetNumLeaf( std::string treeFile );


////// 領域取得関係 //////

  /** 指定面における自リーフの隣接リーフ番号を取得
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(CARTのとき1)
   *  @return 指定面における自リーフの隣接リーフ番号整数配列のポインタ
   */
  virtual
  const int* GetNeighborLeafList( cpm_FaceFlag face, int &num ) const;

  /** 指定面における自リーフの周期境界の隣接リーフ番号を取得
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(CARTのとき1)
   *  @return 指定面における自リーフの周期境界の隣接リーフ番号整数配列のポインタ
   */
  virtual
  const int* GetPeriodicLeafList( cpm_FaceFlag face, int &num ) const;

  /** 指定面における自リーフの隣接ランク番号を取得
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(CARTのとき1)
   *  @return 指定面における自リーフの隣接ランク番号整数配列のポインタ
   */
  virtual
  const int* GetNeighborRankList( cpm_FaceFlag face, int &num ) const;

  /** 指定面における自リーフの周期境界の隣接ランク番号を取得
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(CARTのとき1)
   *  @return 指定面における自リーフの周期境界の隣接ランク番号整数配列のポインタ
   */
  virtual
  const int* GetPeriodicRankList( cpm_FaceFlag face, int &num ) const;

  /** 指定面におけるレベル差を取得
   *  @param[in]  face 面方向
   *  @return     レベル差(0:同じレベル, 1:fine, -1:coarse)
   */
  virtual
  int GetNeighborLevelDiff( cpm_FaceFlag face ) const;

  /** 自リーフの境界が外部境界かどうかを判定
   *  @param[in] face  面方向
   *  @retval    true  外部境界
   *  @retval    false 外部境界でない
   */
  virtual
  bool IsOuterBoundary( cpm_FaceFlag face ) const;

  /** 自リーフの境界が内部境界(隣が不活性ドメイン)かどうかを判定
   *  @param[in] face  面方向
   *  @retval    true  内部境界
   *  @retval    false 内部境界でない
   */
  virtual
  bool IsInnerBoundary( cpm_FaceFlag face ) const;



////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


protected:
  /**** 木情報 ****/
  BCMFileIO::OctHeader m_octHeader;  ///< 木情報ファイルのヘッダー情報
  BCMOctree *m_octree;               ///< 生成された木情報
  Node      *m_node;                 ///< 自ランクが担当するリーフノード
  int        m_leafID;               ///< リーフID

  /**** 隣接情報 ****/
  const NeighborInfo *m_neighborInfo; ///< BCMOctreeから生成した隣接情報
  int m_neighborLeafID_LMR[6][4]; ///< 隣接リーフ番号(外部境界は負の値)
  int m_periodicLeafID_LMR[6][4]; ///< 周期境界の隣接リーフ番号
  int m_neighborRankID_LMR[6][4]; ///< 隣接ランク番号(外部境界は負の値)
  int m_periodicRankID_LMR[6][4]; ///< 周期境界の隣接ランク番号
  int m_neighborLevelDiff[6]; ///< 隣接リーフとのレベル差(-1/0/1)

};

#endif /* _CPM_VOXELINFO_LMR_H_ */
