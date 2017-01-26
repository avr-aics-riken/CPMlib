/*
 * CPMlib - Cartesian Partition Manager Library
 *
 * Copyright (C) 2012-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2014-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   cpm_ParaManagerLMR.cpp
 * パラレルマネージャクラスのソースファイル
 * @date   2012/05/31
 */

#include "cpm_ParaManagerLMR.h"

////////////////////////////////////////////////////////////////////////////////
// インスタンスの取得
cpm_ParaManagerLMR*
cpm_ParaManagerLMR::get_instance()
{
  // 宣言
  static cpm_ParaManagerLMR instance;

  // ポインタ
  return &instance;
}

////////////////////////////////////////////////////////////////////////////////
// インスタンスの取得
cpm_ParaManagerLMR*
cpm_ParaManagerLMR::get_instance(int &argc, char**& argv)
{
#if 0
  // 宣言
  static cpm_ParaManagerLMR instance;

  // MPI_Init
  if( instance.Initialize(argc,argv) != CPM_SUCCESS )
  {
    return NULL;
  }

  // ポインタ
  return &instance;
#else
  // ポインタ取得
  cpm_ParaManagerLMR *instance = cpm_ParaManagerLMR::get_instance();

  // MPI_Init
  if( instance->Initialize(argc,argv) != CPM_SUCCESS )
  {
    return NULL;
  }

  // ポインタ
  return instance;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_ParaManagerLMR::cpm_ParaManagerLMR()
  : cpm_BaseParaManager()
{
  // 領域管理情報マップの削除、クリア
  {
    VoxelInfoMapLMR::iterator it  = m_voxelInfoMap.begin();
    VoxelInfoMapLMR::iterator ite = m_voxelInfoMap.end();
    for( ; it!=ite; it++ )
    {
      LeafMap::iterator ls = it->second.begin();
      LeafMap::iterator le = it->second.end();
      for( ; ls!=le; ls++ )
      {
        if( ls->second ) delete ls->second;
      }
    }
    m_voxelInfoMap.clear();
  }

  // 領域分割タイプ
  m_domainType = CPM_DOMAIN_LMR;
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_ParaManagerLMR::~cpm_ParaManagerLMR()
{
  // 袖通信バッファ情報の削除、クリア
  BndCommInfoMap* pBndCommInfoMapList[6] = { &m_bndCommInfoMapMX
                                           , &m_bndCommInfoMapPX
                                           , &m_bndCommInfoMapMY
                                           , &m_bndCommInfoMapPY
                                           , &m_bndCommInfoMapMZ
                                           , &m_bndCommInfoMapPZ
                                           };
  for( int i=0;i<6;i++ )
  {
    BndCommInfoMap* pBndCommInfoMap = pBndCommInfoMapList[i];
    for( BndCommInfoMap::iterator itP=pBndCommInfoMap->begin();itP!=pBndCommInfoMap->end();itP++ )
    {
      LeafCommInfoMap &CommInfoMap = itP->second;
      for( LeafCommInfoMap::iterator it=CommInfoMap.begin();it!=CommInfoMap.end();it++ )
      {
        delete it->second;
      }
    }
  }
}

  
////////////////////////////////////////////////////////////////////////////////
// LMR用の領域分割
cpm_ErrorCode
cpm_ParaManagerLMR::VoxelInit_LMR( std::string treeFile
                              , size_t maxVC, size_t maxN, int procGrpNo )
{
  cpm_ErrorCode ret = CPM_SUCCESS;
  cpm_VoxelInfoLMR *voxelInfo = NULL;

  // 既に領域分割済みか
  if( m_voxelInfoMap.find(procGrpNo) != m_voxelInfoMap.end() )
  {
    return CPM_ERROR_ALREADY_VOXELINIIT;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm( procGrpNo );
  if( IsCommNull( comm ) )
  {
    return CPM_ERROR_MPI_INVALID_COMM;
  }

  // 定義点管理マップに登録
  DefPointMap::iterator it = m_defPointMap.find(procGrpNo);
  if( it == m_defPointMap.end() ) {
    if( !m_defPointMap.insert(std::make_pair(procGrpNo, CPM_DEFPOINTTYPE_FVM)).second ) {
      return CPM_ERROR_INSERT_DEFPOINTTYPEMAP;
    }
  }

  // 領域分割情報の生成
  LeafMap leafMap;
  if( (ret = cpm_VoxelInfoLMR::Init( comm, treeFile, leafMap )) != CPM_SUCCESS )
  {
    Abort(ret);
    return ret;
  }

  // VOXEL空間マップに登録
  if( !m_voxelInfoMap.insert(std::make_pair(procGrpNo, leafMap)).second )
  {
    for( LeafMap::iterator it=leafMap.begin();it!=leafMap.end();it++ )
    {
      if( it->second )
      {
        delete it->second;
      }
    }
    m_procGrpList.pop_back();
    return CPM_ERROR_INSERT_VOXELMAP;
  }

  // LMR用の袖通信情報を生成
  SetBndCommBuffer( maxVC, maxN, procGrpNo );

  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// 木情報ファイルからリーフ数を取得する
int
cpm_ParaManagerLMR::GetNumLeaf( std::string treeFile )
{
   return cpm_VoxelInfoLMR::GetNumLeaf( treeFile );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファサイズの取得
size_t
cpm_ParaManagerLMR::GetBndCommBufferSize( int procGrpNo )
{
  size_t mem = 0;

  BndCommInfoMap* pBndCommInfoMapList[6] = { &m_bndCommInfoMapMX
                                           , &m_bndCommInfoMapPX
                                           , &m_bndCommInfoMapMY
                                           , &m_bndCommInfoMapPY
                                           , &m_bndCommInfoMapMZ
                                           , &m_bndCommInfoMapPZ
                                           };
  if( procGrpNo < 0 )
  {
    //全体
    for( int i=0;i<6;i++ )
    {
      BndCommInfoMap* pBndCommInfoMap = pBndCommInfoMapList[i];
      for( BndCommInfoMap::iterator itP=pBndCommInfoMap->begin();itP!=pBndCommInfoMap->end();itP++ )
      {
        LeafCommInfoMap &CommInfoMap = itP->second;

        for( LeafCommInfoMap::iterator it=CommInfoMap.begin();it!=CommInfoMap.end();it++ )
        {
          cpm_LeafCommInfo* pLeafCommInfo = it->second;
          mem += pLeafCommInfo->GetBndCommBufferSize();
        }
      }
    }
  }
  else
  {
    // 指定プロセスグループを検索
    for( int i=0;i<6;i++ )
    {
      BndCommInfoMap* pBndCommInfoMap = pBndCommInfoMapList[i];
      BndCommInfoMap::iterator itP = pBndCommInfoMap->find(procGrpNo);
      LeafCommInfoMap &CommInfoMap = itP->second;

      for( LeafCommInfoMap::iterator it=CommInfoMap.begin();it!=CommInfoMap.end();it++ )
      {
        cpm_LeafCommInfo* pLeafCommInfo = it->second;
        mem += pLeafCommInfo->GetBndCommBufferSize();
      }
    }
  }

  mem *= sizeof(REAL_BUF_TYPE);
  return mem;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファのセット
// LMR用の袖通信情報も生成する
cpm_ErrorCode
cpm_ParaManagerLMR::SetBndCommBuffer( size_t maxVC, size_t maxN, int procGrpNo )
{
  if( maxVC==0 || maxN==0 )
  {
    return CPM_ERROR_BNDCOMM;
  }

  // local voxel size
  const int *sz = GetLocalVoxelSize(procGrpNo);
  if( !sz ) return CPM_ERROR_BNDCOMM_VOXELSIZE;

  // 隣接面方向ごとに情報を生成
  const cpm_FaceFlag face[6] = {X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS};
  BndCommInfoMap* pBndCommInfoMapList[6] = { &m_bndCommInfoMapMX
                                           , &m_bndCommInfoMapPX
                                           , &m_bndCommInfoMapMY
                                           , &m_bndCommInfoMapPY
                                           , &m_bndCommInfoMapMZ
                                           , &m_bndCommInfoMapPZ
                                           };
  size_t sz_face[6][2] = { {sz[1], sz[2]}, {sz[1], sz[2]}
                         , {sz[2], sz[0]}, {sz[2], sz[0]}
                         , {sz[0], sz[1]}, {sz[0], sz[1]} };
  for( int i=0;i<6;i++ )
  {
    BndCommInfoMap *pBndCommInfoMap = pBndCommInfoMapList[i];

    // プロセスグループを検索
    BndCommInfoMap::iterator itP = pBndCommInfoMap->find(procGrpNo);
    if( itP != pBndCommInfoMap->end() )
    {
      // すでに生成されている
      continue;
    }

    // インスタンス
    LeafCommInfoMap commInfoMap;

    // 自ランクに含まれるリーフについてすべて検索する
    std::vector<int> leafIDs = GetLocalLeafIDs(procGrpNo);
    for( int j=0;j<leafIDs.size();j++ )
    {
      // 周期境界フラグ
      bool bPeriodic = false;

      // VoxelInfo
      int leafID = leafIDs[j];
      const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo_byID(leafID, procGrpNo);

      // 隣接リーフとランクを取得
      int num = 0;
      const int *distLeafList = pVoxelInfo->GetNeighborLeafList(face[i], num);
      const int *distRankList = pVoxelInfo->GetNeighborRankList(face[i], num);
      if( num == 0 )
      {
        // 内部の隣接リーフが無いとき、周期境界方向から取得
        distLeafList = pVoxelInfo->GetPeriodicLeafList(face[i], num);
        distRankList = pVoxelInfo->GetPeriodicRankList(face[i], num);
        bPeriodic = true;
      }
      if( num == 0 )
      {
        // 隣接リーフが無い
        continue;
      }

      // レベル差を取得
      int levelDiff = pVoxelInfo->GetNeighborLevelDiff(face[i]);

      // 各面ごとにマップに登録
      for( int n=0;n<num;n++ )
      {
        // 隣接リーフとランク
        int distLeafID = distLeafList[n];
        int distRankID = distRankList[n];

        // stCommInfo
        cpm_LeafCommInfo::stCommInfo *commInfo = new cpm_LeafCommInfo::stCommInfo();
        commInfo->iOwnLeafID  = leafID;
        commInfo->iDistLeafID = distLeafID;
        commInfo->iLevelDiff  = levelDiff;
        commInfo->iFaceIdx    = n;
        commInfo->bPeriodic   = bPeriodic;

        // LeafCommInfoに格納
        cpm_LeafCommInfo *pLeafCommInfo = NULL;
        LeafCommInfoMap::iterator itR = commInfoMap.find(distRankID);
        if( itR == commInfoMap.end() )
        {
          pLeafCommInfo = new cpm_LeafCommInfo(distRankID);
          commInfoMap.insert(LeafCommInfoMap::value_type(distRankID, pLeafCommInfo));
        }
        else
        {
          pLeafCommInfo = itR->second;
        }
        pLeafCommInfo->AddCommInfo(commInfo);
      } //for face(1-4)

    } //for leafIDs

    // マップに登録
    pBndCommInfoMap->insert(BndCommInfoMap::value_type(procGrpNo, commInfoMap));

  } //for 6 (face dir)

  // 通信情報をLeafIDで整列し、送受信/コピー用のバッファを確保する
  // 自身のランク番号 <= 相手のランク番号のとき  自身のリーフ番号で昇順
  // 自身のランク番号 >  相手のランク番号のとき  相手のリーフ番号で昇順
  for( int i=0;i<6;i++ )
  {
    // LeafCommInfoMapの取得
    BndCommInfoMap *pBndCommInfoMap = pBndCommInfoMapList[i];
    BndCommInfoMap::iterator itP = pBndCommInfoMap->find(procGrpNo);
    LeafCommInfoMap &CommInfoMap = itP->second;

    // 通信相手ランクごとにソートする
    for( LeafCommInfoMap::iterator it=CommInfoMap.begin();it!=CommInfoMap.end();it++ )
    {
      int distRank = it->first;
      cpm_LeafCommInfo* pLeafCommInfo = it->second;
      if( m_rankNo <= distRank )
      {
        pLeafCommInfo->Sort(0);
      }
      else
      {
        pLeafCommInfo->Sort(1);
      }

      // 送受信/コピー用のバッファの確保
      pLeafCommInfo->SetBndCommBuffer( m_rankNo, sz_face[i], maxVC, maxN );
    }
  }

#if 1
fflush(stdout);
Barrier(procGrpNo);
{
  char fname[256];
  sprintf( fname, "LMR_CommInfo_rank%d.log", m_rankNo);
  FILE *fp = fopen(fname, "wt");
  const char* faceName[] = {"X_MINUS", "X_PLUS", "Y_MINUS", "Y_PLUS", "Z_MINUS", "Z_PLUS"};
  for( int i=0;i<6;i++ )
  {
    fprintf(fp, "****** face %s CommInfo *************************************************\n", faceName[i]);
    BndCommInfoMap *pBndCommInfoMap = pBndCommInfoMapList[i];
    BndCommInfoMap::iterator itP = pBndCommInfoMap->find(procGrpNo);
    LeafCommInfoMap &CommInfoMap = itP->second;
    for( LeafCommInfoMap::iterator it=CommInfoMap.begin();it!=CommInfoMap.end();it++ )
    {
      int distRank = it->first;
      cpm_LeafCommInfo* pLeafCommInfo = it->second;
      fprintf(fp, "  *distRank = %d\n", distRank);
      for( int j=0;j<pLeafCommInfo->m_vecCommInfo.size();j++ )
      {
        cpm_LeafCommInfo::stCommInfo* commInfo = pLeafCommInfo->m_vecCommInfo[j];
        fprintf(fp, "    *ownLeaf=%d, distLeaf=%d\n", commInfo->iOwnLeafID, commInfo->iDistLeafID);
        fprintf(fp, "     iLevelDiff=%d\n", commInfo->iLevelDiff);
        fprintf(fp, "     iFaceIdx  =%d\n", commInfo->iFaceIdx);
        if( commInfo->bPeriodic )
          fprintf(fp, "     periodic  =true\n");
        else
          fprintf(fp, "     periodic  =false\n");
      }
    }
  }
  fclose(fp);
}
fflush(stdout);
Barrier(procGrpNo);
#endif

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// VOXEL空間マップを検索
const cpm_VoxelInfo*
cpm_ParaManagerLMR::FindVoxelInfo( int procGrpNo )
{
  return FindLeafVoxelInfo(0, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// VOXEL空間マップを検索
const cpm_VoxelInfoLMR*
cpm_ParaManagerLMR::FindLeafVoxelInfo( int leafIndex, int procGrpNo )
{
  std::vector<int> leafIDs = GetLocalLeafIDs(procGrpNo);
  if( leafIDs.size() == 0 )
  {
    return NULL;
  }
  return FindLeafVoxelInfo_byID(leafIDs[leafIndex], procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// VOXEL空間マップを検索
const cpm_VoxelInfoLMR*
cpm_ParaManagerLMR::FindLeafVoxelInfo_byID( int leafID, int procGrpNo )
{
  VoxelInfoMapLMR::iterator it = m_voxelInfoMap.find(procGrpNo);
  if( it == m_voxelInfoMap.end() ) return NULL;

  LeafMap::iterator itLeaf = it->second.find(leafID);
  if( itLeaf == it->second.end() ) return NULL;

  return itLeaf->second;
}

////////////////////////////////////////////////////////////////////////////////
// 全リーフ数を取得する
int
cpm_ParaManagerLMR::GetNumLeaf( int procGrpNo )
{
  int nbuf = GetLocalNumLeaf(procGrpNo);
  int nLeaf = nbuf;
  Allreduce(&nbuf, &nLeaf, 1, MPI_SUM, procGrpNo);
  return nLeaf;
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクが担当するリーフ数を取得する
int
cpm_ParaManagerLMR::GetLocalNumLeaf( int procGrpNo )
{
  std::vector<int> leafIDs = GetLocalLeafIDs(procGrpNo);
  return leafIDs.size();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクが担当するリーフIDリストを取得する
std::vector<int>
cpm_ParaManagerLMR::GetLocalLeafIDs( int procGrpNo )
{
  std::vector<int> leafIDs;
  VoxelInfoMapLMR::iterator it = m_voxelInfoMap.find(procGrpNo);
  if( it == m_voxelInfoMap.end() )
  {
    return leafIDs;
  }

  for( LeafMap::iterator itL=it->second.begin();itL!=it->second.end();itL++ )
  {
    leafIDs.push_back(itL->second->m_leafID);
  }

  return leafIDs;
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフのリーフIDを取得する
int
cpm_ParaManagerLMR::GetLeafID( int leafIndex, int procGrpNo )
{
  std::vector<int> leafIDs = GetLocalLeafIDs(procGrpNo);
  return leafIDs[leafIndex];
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクが担当するリーフIDのインデクスを取得する
int
cpm_ParaManagerLMR::GetLocalLeafIndex_byID( int leafID, int procGrpNo )
{
  std::vector<int> leafIDs = GetLocalLeafIDs(procGrpNo);
  for( size_t i=0;i<leafIDs.size();i++ )
  {
    if( leafIDs[i] == leafID )
    {
      return (int)i;
    }
  }
  return -1;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割数を取得
const int*
cpm_ParaManagerLMR::GetDivNum( int leafIndex, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetDivNum();
}

////////////////////////////////////////////////////////////////////////////////
// ピッチを取得
const double*
cpm_ParaManagerLMR::GetPitch( int leafIndex, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPitch();
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの空間原点を取得
const double*
cpm_ParaManagerLMR::GetLocalOrigin( int leafIndex, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの空間サイズを取得
const double*
cpm_ParaManagerLMR::GetLocalRegion( int leafIndex, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの領域分割位置を取得
const int*
cpm_ParaManagerLMR::GetDivPos( int leafIndex, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetDivPos();
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの始点VOXELの全体空間でのインデクスを取得
const int*
cpm_ParaManagerLMR::GetVoxelHeadIndex( int leafIndex, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetVoxelHeadIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの始点頂点の全体空間でのインデクスを取得
const int*
cpm_ParaManagerLMR::GetNodeHeadIndex( int leafIndex, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNodeHeadIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの始点VOXELまたは頂点の全体空間でのインデクスを取得
const int*
cpm_ParaManagerLMR::GetArrayHeadIndex( int leafIndex, int procGrpNo )
{
  //定義点がVOXELのときVOXELのインデックスを取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FVM )
  {
    return GetVoxelHeadIndex(leafIndex, procGrpNo);
  }

  //定義点がNODEのとき頂点のインデックスを取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FDM )
  {
    return GetNodeHeadIndex(leafIndex, procGrpNo);
  }

  return NULL;
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの終点VOXELの全体空間でのインデクスを取得
const int*
cpm_ParaManagerLMR::GetVoxelTailIndex( int leafIndex, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetVoxelTailIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの終点頂点の全体空間でのインデクスを取得
const int*
cpm_ParaManagerLMR::GetNodeTailIndex( int leafIndex, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNodeTailIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの終点VOXELまたは頂点の全体空間でのインデクスを取得
const int*
cpm_ParaManagerLMR::GetArrayTailIndex( int leafIndex, int procGrpNo )
{
  //定義点がVOXELのときVOXELのインデックスを取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FVM )
  {
    return GetVoxelTailIndex(leafIndex, procGrpNo);
  }

  //定義点がNODEのとき頂点のインデックスを取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FDM )
  {
    return GetNodeTailIndex(leafIndex, procGrpNo);
  }

  return NULL;
}

////////////////////////////////////////////////////////////////////////////////
// 指定面における自リーフの隣接リーフ番号を取得
const int*
cpm_ParaManagerLMR::GetNeighborLeafList( int leafIndex, cpm_FaceFlag face, int &num, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNeighborLeafList( face, num );
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの指定面における自リーフの周期境界の隣接リーフ番号を取得
const int*
cpm_ParaManagerLMR::GetPeriodicLeafList( int leafIndex, cpm_FaceFlag face, int &num, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPeriodicLeafList( face, num );
}

////////////////////////////////////////////////////////////////////////////////
// 指定面における自リーフの隣接ランク番号を取得
const int*
cpm_ParaManagerLMR::GetNeighborRankList( int leafIndex, cpm_FaceFlag face, int &num, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNeighborRankList( face, num );
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの指定面における自リーフの周期境界の隣接ランク番号を取得
const int*
cpm_ParaManagerLMR::GetPeriodicRankList( int leafIndex, cpm_FaceFlag face, int &num, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPeriodicRankList( face, num );
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの指定面におけるレベル差を取得
int
cpm_ParaManagerLMR::GetNeighborLevelDiff( int leafIndex, cpm_FaceFlag face, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNeighborLevelDiff( face );
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの境界が外部境界かどうかを判定
bool
cpm_ParaManagerLMR::IsOuterBoundary( int leafIndex, cpm_FaceFlag face, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->IsOuterBoundary( face );
}

////////////////////////////////////////////////////////////////////////////////
// 指定リーフの境界が内部境界(隣が不活性ドメイン)かどうかを判定
bool
cpm_ParaManagerLMR::IsInnerBoundary( int leafIndex, cpm_FaceFlag face, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfoLMR *pVoxelInfo = FindLeafVoxelInfo( leafIndex, procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->IsInnerBoundary( face );
}

