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
 * @file   cpm_ParaManager.cpp
 * カーテシアン用パラレルマネージャクラスのソースファイル
 * @date   2015/03/27
 */

#include "cpm_ParaManager.h"
#include "cpm_VoxelInfoCART.h"

////////////////////////////////////////////////////////////////////////////////
// 唯一のインスタンスの取得
cpm_ParaManager*
cpm_ParaManager::get_instance()
{
  // 宣言
  static cpm_ParaManager instance;

  // ポインタ
  return &instance;
}

////////////////////////////////////////////////////////////////////////////////
// 唯一のインスタンスの取得(initialize処理も実行)
cpm_ParaManager*
cpm_ParaManager::get_instance(int &argc, char**& argv)
{
#if 0
  // 宣言
  static cpm_ParaManager instance;

  // MPI_Init
  if( instance.Initialize(argc,argv) != CPM_SUCCESS )
  {
    return NULL;
  }

  // ポインタ
  return &instance;
#else
  // ポインタ取得
  cpm_ParaManager *instance = cpm_ParaManager::get_instance();

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
cpm_ParaManager::cpm_ParaManager()
  : cpm_BaseParaManager()
{
  // 領域分割タイプ
  m_domainType = CPM_DOMAIN_CARTESIAN;

  // 領域管理情報マップのクリア
  m_voxelInfoMap.clear();

  // 袖通信バッファ情報のクリア
  m_bndCommInfoMap.clear();
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_ParaManager::~cpm_ParaManager()
{

  // 領域管理情報マップの削除、クリア
  {
    VoxelInfoMap::iterator it  = m_voxelInfoMap.begin();
    VoxelInfoMap::iterator ite = m_voxelInfoMap.end();
    for( ; it!=ite; it++ )
    {
      if( it->second ) delete it->second;
    }
    m_voxelInfoMap.clear();
  }

  // 袖通信バッファ情報の削除、クリア
  {
    BndCommInfoMap::iterator it  = m_bndCommInfoMap.begin();
    BndCommInfoMap::iterator ite = m_bndCommInfoMap.end();
    for( ; it!=ite; it++ )
    {
      if( it->second ) delete it->second;
    }
    m_bndCommInfoMap.clear();
  }
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(FVM用)
cpm_ErrorCode
cpm_ParaManager::VoxelInit( cpm_GlobalDomainInfo* domainInfo, size_t maxVC, size_t maxN, int procGrpNo )
{
  cpm_ErrorCode ret;
  cpm_VoxelInfoCART *voxelInfo = NULL;

  // 入力引数のチェック
  if( !domainInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

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

  // commの並列数を取得
  int nRank = 0;
  MPI_Comm_size( comm, &nRank );

  // 領域情報のチェック
  // 活性サブドメイン配列が空のとき、全領域が活性サブドメインになる
  // その場合、CheckData関数内で活性サブドメイン情報を生成する
  if( (ret = domainInfo->CheckData( nRank )) != CPM_SUCCESS )
  {
    return ret;
  }

  // プロセス数と活性サブドメイン数のチェック
  if( nRank != domainInfo->GetSubdomainNum() )
  {
    Abort(CPM_ERROR_MISMATCH_NP_SUBDOMAIN);
    return CPM_ERROR_MISMATCH_NP_SUBDOMAIN;
  }

  // インスタンス
  voxelInfo = new cpm_VoxelInfoCART();
  if( !voxelInfo )
  {
    Abort(CPM_ERROR_INVALID_PTR);
    return CPM_ERROR_INVALID_PTR;
  }

  // 領域分割情報の生成
  if( (ret = voxelInfo->Init( comm, domainInfo )) != CPM_SUCCESS )
  {
    delete voxelInfo;
    Abort(ret);
    return ret;
  }

  // VOXEL空間マップに登録
  if( !m_voxelInfoMap.insert(std::make_pair(procGrpNo, voxelInfo)).second )
  {
    m_procGrpList.pop_back();
    delete voxelInfo;
    return CPM_ERROR_INSERT_VOXELMAP;
  }

  // 袖通信バッファの設定
  SetBndCommBuffer( maxVC, maxN, procGrpNo );

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(領域分割数を指定)(FVM用)
cpm_ErrorCode
cpm_ParaManager::VoxelInit( int div[3], int vox[3], double origin[3], double region[3]
                              , size_t maxVC, size_t maxN, cpm_DivPolicy divPolicy, int procGrpNo )
{
  // 入力値のチェック
  if( vox[0] <= 0 || vox[1] <= 0 || vox[2] <= 0 )
  {
    return CPM_ERROR_INVALID_VOXELSIZE;
  }
  if( region[0] <= 0.0 || region[1] <= 0.0 || region[2] <= 0.0 )
  {
    return CPM_ERROR_INVALID_REGION;
  }

  // 領域分割数
  if( div[0]<=0 || div[1]<=0 || div[2]<=0 )
  {
    // ランク数
    int nrank = GetNumRank( procGrpNo );

    // 領域分割数の決定
    div[0] = div[1] = div[2] = 0;
    cpm_ErrorCode ret;
    if( divPolicy == DIV_COMM_SIZE )
    {
      ret = DecideDivPattern_CommSize( nrank, vox, div );
    }
    else
    {
      ret = DecideDivPattern_Cube( nrank, vox, div );
    }
    if( ret != CPM_SUCCESS )
    {
      return ret;
    }
  }

  //ピッチを計算
  double pitch[3];
  for( int i=0;i<3;i++ )
  {
    pitch[i] = region[i] / double(vox[i]);
  }

  // DomainInfoを生成
  cpm_GlobalDomainInfo dInfo;
  dInfo.SetOrigin( origin );
  dInfo.SetPitch ( pitch );
  dInfo.SetRegion( region );
  dInfo.SetVoxNum( vox );
  dInfo.SetDivNum( div );

  // 共通の処理
  return VoxelInit( &dInfo, maxVC, maxN, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(プロセスグループのランク数で自動領域分割)(FVM用)
cpm_ErrorCode
cpm_ParaManager::VoxelInit( int vox[3], double origin[3], double region[3]
                              , size_t maxVC, size_t maxN, cpm_DivPolicy divPolicy
                              , int procGrpNo )
{
  // 領域分割
  int div[3] = {0,0,0};
  return VoxelInit( div, vox, origin, region, maxVC, maxN, divPolicy, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(領域分割数を指定)(FDM用)
cpm_ErrorCode
cpm_ParaManager::NodeInit ( int div[3], int nod[3], double origin[3], double region[3]
                              , size_t maxVC, size_t maxN, cpm_DivPolicy divPolicy, int procGrpNo )
{

  cpm_ErrorCode ret;

  // 定義点管理マップに登録
  if( !m_defPointMap.insert(std::make_pair(procGrpNo, CPM_DEFPOINTTYPE_FDM)).second )
  {
    return CPM_ERROR_INSERT_DEFPOINTTYPEMAP;
  }

  // VOXEL数のセット
  int vox[3];
  vox[0] = nod[0] - 1;
  vox[1] = nod[1] - 1;
  vox[2] = nod[2] - 1;
  ret = VoxelInit( div, vox, origin, region, maxVC, maxN, divPolicy, procGrpNo );
  if( ret != CPM_SUCCESS ) {
    if(      ret == CPM_ERROR_ALREADY_VOXELINIIT ) {
      return CPM_ERROR_ALREADY_NODEINIT;
    }
    else if( ret == CPM_ERROR_INVALID_VOXELSIZE  )  {
      return CPM_ERROR_INVALID_NODESIZE;
    }
    return ret;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(プロセスグループのランク数で自動領域分割)(FDM用)
cpm_ErrorCode
cpm_ParaManager::NodeInit ( int nod[3], double origin[3], double region[3]
                              , size_t maxVC, size_t maxN, cpm_DivPolicy divPolicy
                              , int procGrpNo )
{

  // 領域分割
  int div[3] = {0, 0, 0};
  return NodeInit( div, nod, origin, region, maxVC, maxN, divPolicy, procGrpNo );

}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(ActiveSubdomainファイル、領域分割数を指定)(FVM用)
cpm_ErrorCode
cpm_ParaManager::VoxelInit_Subdomain( int div[3], int vox[3], double origin[3], double region[3]
                                        , std::string subDomainFile
                                        , size_t maxVC, size_t maxN, int procGrpNo )
{
  // 入力値のチェック
  if( vox[0] <= 0 || vox[1] <= 0 || vox[2] <= 0 )
  {
    return CPM_ERROR_INVALID_VOXELSIZE;
  }
  if( region[0] <= 0.0 || region[1] <= 0.0 || region[2] <= 0.0 )
  {
    return CPM_ERROR_INVALID_REGION;
  }

  // DomainInfoを生成
  cpm_GlobalDomainInfo dInfo;

  // ActiveSubdomainファイルを読み込み
  int divSubdomain[3] = {0,0,0};
  std::vector<cpm_ActiveSubdomainInfo> subDomainInfo;
  cpm_ErrorCode ret = cpm_GlobalDomainInfo::ReadActiveSubdomainFile
               (         subDomainFile, subDomainInfo, divSubdomain );
  if( ret != CPM_SUCCESS )
  {
    return ret;
  }

  // ActiveSubdomain数
  int ndom = subDomainInfo.size();

  // ランク数
  int nrank = GetNumRank( procGrpNo );
  if( nrank != ndom )
  {
    return CPM_ERROR_MISMATCH_NP_SUBDOMAIN;
  }

  // 領域分割数
  if( div[0]<=0 || div[1]<=0 || div[2]<=0 )
  {
    // 領域分割数=ActiveSubdomainファイルの領域分割数
    div[0] = divSubdomain[0];
    div[1] = divSubdomain[1];
    div[2] = divSubdomain[2];
  }
  else
  {
    if( div[0] != divSubdomain[0] ||
        div[1] != divSubdomain[1] ||
        div[2] != divSubdomain[2] )
    {
      return CPM_ERROR_MISMATCH_DIV_SUBDOMAIN;
    }
  }

  //ピッチを計算
  double pitch[3];
  for( int i=0;i<3;i++ )
  {
    pitch[i] = region[i] / double(vox[i]);
  }

  // DomainInfoに値をセット
  dInfo.SetOrigin( origin );
  dInfo.SetPitch ( pitch );
  dInfo.SetRegion( region );
  dInfo.SetVoxNum( vox );
  dInfo.SetDivNum( div );
  for( size_t i=0;i<subDomainInfo.size();i++ )
  {
    dInfo.AddSubdomain(subDomainInfo[i]);
  }

  // 共通の処理
  return VoxelInit( &dInfo, maxVC, maxN, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(ActiveSubdomainファイルを指定)(FVM用)
cpm_ErrorCode
cpm_ParaManager::VoxelInit_Subdomain( int vox[3], double origin[3], double region[3]
                                        , std::string subDomainFile
                                        , size_t maxVC, size_t maxN, int procGrpNo )
{
  // 領域分割
  int div[3] = {0,0,0};
  return VoxelInit_Subdomain( div, vox, origin, region, subDomainFile
                            , maxVC, maxN, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(ActiveSubdomainファイル、領域分割数を指定)(FDM用)
cpm_ErrorCode
cpm_ParaManager::NodeInit_Subdomain( int div[3], int nod[3], double origin[3], double region[3]
                                        , std::string subDomainFile
                                        , size_t maxVC, size_t maxN, int procGrpNo )
{

  cpm_ErrorCode ret;

  // 定義点管理マップに登録
  if( !m_defPointMap.insert(std::make_pair(procGrpNo, CPM_DEFPOINTTYPE_FDM)).second )
  {
    return CPM_ERROR_INSERT_DEFPOINTTYPEMAP;
  }

  // VOXEL数のセット
  int vox[3];
  vox[0] = nod[0] - 1;
  vox[1] = nod[1] - 1;
  vox[2] = nod[2] - 1;
  ret = VoxelInit_Subdomain( div, vox, origin, region, subDomainFile
                           , maxVC, maxN, procGrpNo );
  if( ret != CPM_SUCCESS ) {
    return ret;
  }

  return CPM_SUCCESS;

}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(ActiveSubdomainファイルを指定)(FDM用)
cpm_ErrorCode
cpm_ParaManager::NodeInit_Subdomain( int nod[3], double origin[3], double region[3]
                                        , std::string subDomainFile
                                        , size_t maxVC, size_t maxN, int procGrpNo )
{
  // 領域分割
  int div[3] = {0,0,0};
  return NodeInit_Subdomain( div, nod, origin, region, subDomainFile
                            , maxVC, maxN, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 並列プロセス数からI,J,K方向の分割数を取得する
// 通信面のトータルサイズが小さい分割パターンを採用する
cpm_ErrorCode
cpm_ParaManager::DecideDivPattern_CommSize( int divNum
                                              , int voxSize[3]
                                              , int divPttn[3]
                                              ) const
{
  if( !voxSize || !divPttn )
  {
    return CPM_ERROR_INVALID_PTR;
  }
  if( (voxSize[0]==0) || (voxSize[1]==0) || (voxSize[2]==0) )
  {
    return CPM_ERROR_INVALID_VOXELSIZE;
  }
  if( divNum <= 1 ){
    divPttn[0] = divPttn[1] = divPttn[2] = 1;
    return CPM_SUCCESS;
  }

  divPttn[0] = divPttn[1] = divPttn[2] = 0;

  unsigned long long minCommSize = 0;

  unsigned long long divNumll = divNum;
  unsigned long long voxSizell[3] = {0, 0, 0};
  unsigned long long divPttnll[3] = {0, 0, 0};
  voxSizell[0] = voxSize[0];
  voxSizell[1] = voxSize[1];
  voxSizell[2] = voxSize[2];

  bool flag = false;
  unsigned long long i, j, k;
  for(i=1; i<=divNumll; i++)
  {
    if( divNumll%i != 0 ) continue;
    if( voxSizell[0] < i ) break;
    unsigned long long jmax = divNumll/i;
    for(j=1; j<=jmax; j++)
    {
      if( jmax%j != 0 ) continue;
      if( voxSizell[1] < j ) break;

      k = jmax/j;
      if( voxSizell[2] < k ) continue;

      unsigned long long commSize;
      if( (commSize=CalcCommSize(i, j, k, voxSizell)) == 0 ) break;

      if( !flag )
      {
        divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
        minCommSize = commSize;
        flag = true;
      }
      else if( commSize < minCommSize )
      {
        divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
        minCommSize = commSize;
      }
    }
  }

  divPttn[0] = divPttnll[0];
  divPttn[1] = divPttnll[1];
  divPttn[2] = divPttnll[2];

  if( (divPttn[0]==0) || (divPttn[1]==0) || (divPttn[2]==0) )
  {
    return CPM_ERROR_DECIDE_DIV_PATTERN;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// I,J,K分割を行った時の通信点数の総数を取得する
unsigned long long
cpm_ParaManager::CalcCommSize( unsigned long long iDiv
                                 , unsigned long long jDiv
                                 , unsigned long long kDiv
                                 , unsigned long long voxSize[3]
                                 ) const
{
  if( (iDiv==0) || (jDiv==0) || (kDiv==0) ) return 0;
  if( !voxSize ) return 0;

  unsigned long long Len[3];
  Len[0] = voxSize[0] / iDiv; if( Len[0] == 0 ) return 0;
  Len[1] = voxSize[1] / jDiv; if( Len[1] == 0 ) return 0;
  Len[2] = voxSize[2] / kDiv; if( Len[2] == 0 ) return 0;

  unsigned long long commFace[3];
  if( iDiv != 1) commFace[0] = Len[1]*Len[2]*(iDiv-1);
  else commFace[0] = 0;
  if( jDiv != 1) commFace[1] = Len[2]*Len[0]*(jDiv-1);
  else commFace[1] = 0;
  if( kDiv != 1) commFace[2] = Len[0]*Len[1]*(kDiv-1);
  else commFace[2] = 0;

  return (commFace[0] + commFace[1] + commFace[2]);
}

////////////////////////////////////////////////////////////////////////////////
// 並列プロセス数からI,J,K方向の分割数を取得する
// １つのサブドメインが立方体に一番近い分割パターンを採用する
cpm_ErrorCode
cpm_ParaManager::DecideDivPattern_Cube( int divNum
                                          , int voxSize[3]
                                          , int divPttn[3]
                                          ) const
{
  if( !voxSize || !divPttn )
  {
    return CPM_ERROR_INVALID_PTR;
  }
  if( (voxSize[0]==0) || (voxSize[1]==0) || (voxSize[2]==0) )
  {
    return CPM_ERROR_INVALID_VOXELSIZE;
  }
  if( divNum <= 1 ){
    divPttn[0] = divPttn[1] = divPttn[2] = 1;
    return CPM_SUCCESS;
  }

  divPttn[0] = divPttn[1] = divPttn[2] = 0;

  unsigned long long minVoxDiff = 0;

  unsigned long long divNumll = divNum;
  unsigned long long voxSizell[3] = {0, 0, 0};
  unsigned long long divPttnll[3] = {0, 0, 0};
  voxSizell[0] = voxSize[0];
  voxSizell[1] = voxSize[1];
  voxSizell[2] = voxSize[2];

  bool flag = false;
  unsigned long long i, j, k;
  for(i=1; i<=divNumll; i++)
  {
    if( divNumll%i != 0 ) continue;
    if( voxSizell[0] < i ) break;
    unsigned long long jmax = divNumll/i;
    for(j=1; j<=jmax; j++)
    {
      if( jmax%j != 0 ) continue;
      if( voxSizell[1] < j ) break;

      k = jmax/j;
      if( voxSizell[2] < k ) continue;

      long long voxDiff;
      if( (voxDiff=CheckCube(i, j, k, voxSizell)) < 0 ) break;

      if( !flag )
      {
        divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
        minVoxDiff = voxDiff;
        flag = true;
      }
      else if( voxDiff < minVoxDiff )
      {
        divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
        minVoxDiff = voxDiff;
      }
    }
  }

  divPttn[0] = divPttnll[0];
  divPttn[1] = divPttnll[1];
  divPttn[2] = divPttnll[2];

  if( (divPttn[0]==0) || (divPttn[1]==0) || (divPttn[2]==0) )
  {
    return CPM_ERROR_DECIDE_DIV_PATTERN;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// I,J,K分割を行った時のI,J,Kボクセル数の最大/最小の差を取得する
long long
cpm_ParaManager::CheckCube( unsigned long long iDiv
                              , unsigned long long jDiv
                              , unsigned long long kDiv
                              , unsigned long long voxSize[3]
                              ) const
{
  if( (iDiv==0) || (jDiv==0) || (kDiv==0) ) return -1;
  if( !voxSize ) return -1;

  unsigned long long Len[3];
  Len[0] = voxSize[0] / iDiv; if( Len[0] == 0 ) return -1;
  Len[1] = voxSize[1] / jDiv; if( Len[1] == 0 ) return -1;
  Len[2] = voxSize[2] / kDiv; if( Len[2] == 0 ) return -1;

  unsigned long long minVox = (Len[0]<Len[1]?Len[0]:Len[1]);
  minVox = (minVox<Len[2]?minVox:Len[2]);
  unsigned long long maxVox = (Len[0]>Len[1]?Len[0]:Len[1]);
  maxVox = (maxVox>Len[2]?maxVox:Len[2]);

  return (maxVox-minVox);
}

////////////////////////////////////////////////////////////////////////////////
// VOXEL空間マップを検索
const cpm_VoxelInfo*
cpm_ParaManager::FindVoxelInfo( int procGrpNo )
{
  VoxelInfoMap::iterator it = m_voxelInfoMap.find(procGrpNo);
  if( it == m_voxelInfoMap.end() ) return NULL;
  return it->second;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割数を取得
const int*
cpm_ParaManager::GetDivNum( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetDivNum();
}

////////////////////////////////////////////////////////////////////////////////
// ピッチを取得
const double*
cpm_ParaManager::GetPitch( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPitch();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの空間原点を取得
const double*
cpm_ParaManager::GetLocalOrigin( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの空間サイズを取得
const double*
cpm_ParaManager::GetLocalRegion( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの領域分割位置を取得
const int*
cpm_ParaManager::GetDivPos( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetDivPos();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの始点VOXELの全体空間でのインデクスを取得
const int*
cpm_ParaManager::GetVoxelHeadIndex( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetVoxelHeadIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの始点頂点の全体空間でのインデクスを取得
const int*
cpm_ParaManager::GetNodeHeadIndex( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNodeHeadIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの始点VOXELまたは頂点の全体空間でのインデクスを取得
// FVMのときはボクセル、FDMのときは頂点での始点インデックスを取得
const int*
cpm_ParaManager::GetArrayHeadIndex( int procGrpNo )
{
  //定義点がVOXELのときVOXELのインデックスを取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FVM ) {
    return GetVoxelHeadIndex(procGrpNo);
  }

  //定義点がNODEのとき頂点のインデックスを取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FDM ) {
    return GetNodeHeadIndex(procGrpNo);
  }

  return NULL;
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの終点VOXELの全体空間でのインデクスを取得
const int*
cpm_ParaManager::GetVoxelTailIndex( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetVoxelTailIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの終点頂点の全体空間でのインデクスを取得
const int*
cpm_ParaManager::GetNodeTailIndex( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNodeTailIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの終点VOXELまたは頂点の全体空間でのインデクスを取得
// FVMのときはボクセル、FDMのときは頂点での終点インデックスを取得
const int*
cpm_ParaManager::GetArrayTailIndex( int procGrpNo )
{
  //定義点がVOXELのときVOXELのインデックスを取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FVM ) {
    return GetVoxelTailIndex(procGrpNo);
  }

  //定義点がNODEのとき頂点のインデックスを取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FDM ) {
    return GetNodeTailIndex(procGrpNo);
  }

  return NULL;
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの隣接ランク番号を取得
const int*
cpm_ParaManager::GetNeighborRankID( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNeighborRankID();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの周期境界の隣接ランク番号を取得
const int*
cpm_ParaManager::GetPeriodicRankID( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPeriodicRankID();
}

////////////////////////////////////////////////////////////////////////////////
// 指定idを含む全体ボクセル空間のインデクス範囲を取得
bool
cpm_ParaManager::GetBndIndexExtGc( int id, int *array, int vc
                                     , int &ista, int &jsta, int &ksta
                                     , int &ilen, int &jlen, int &klen
                                     , int procGrpNo, CPM_PADDING padding )
{
  //ローカルボクセル数

  const int *sz = GetLocalArraySize(procGrpNo);

  if( !sz ) return false;
  int imax = sz[0];
  int jmax = sz[1];
  int kmax = sz[2];

  // パディングサイズ
  int pad_size[3];
  GetPaddingSize(CPM_ARRAY_S3D, sz, vc, pad_size);

  // 共有関数
  return GetBndIndexExtGc( id, array
                         , imax, jmax, kmax, vc
                         , ista, jsta, ksta, ilen, jlen, klen, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 指定idを含む全体ボクセル空間のインデクス範囲を取得
bool
cpm_ParaManager::GetBndIndexExtGc( int id, int *array
                                     , int imax, int jmax, int kmax, int vc
                                     , int &ista, int &jsta, int &ksta
                                     , int &ilen, int &jlen, int &klen
                                     , int pad_size[3], int procGrpNo )
{
  // 全体ボクセルサイズ

  const int *wsz = GetGlobalArraySize(procGrpNo);
  if( !wsz ) return false;

  // 自ランク内のid範囲を取得
  int my_sta[3] = {wsz[0]+vc-1, wsz[1]+vc-1, wsz[2]+vc-1};
  int my_end[3] = {-vc, -vc, -vc};
  for( int k=-vc;k<kmax+vc;k++ ){
  for( int j=-vc;j<jmax+vc;j++ ){
  for( int i=-vc;i<imax+vc;i++ ){
    long long idx = _IDX_S3D(i,j,k,imax,jmax,kmax,vc);
    if( array[idx] == id )
    {
      if( i < my_sta[0] ) my_sta[0] = i;
      if( j < my_sta[1] ) my_sta[1] = j;
      if( k < my_sta[2] ) my_sta[2] = k;
      if( i > my_end[0] ) my_end[0] = i;
      if( j > my_end[1] ) my_end[1] = j;
      if( k > my_end[2] ) my_end[2] = k;
    }
  }}}

  // 取得した自ランク範囲を全体インデクスに変換
  const int *hidx = GetArrayHeadIndex(procGrpNo);
  if( !hidx ) return false;
  for( int i=0;i<3;i++ )
  {
    if( my_sta[i] > my_end[i] ) continue;
    my_sta[i] += hidx[i];
    my_end[i] += hidx[i];
  }

  // 全ランクのmin/maxを取得
  int idxsta[3], idxend[3];
  for( int i=0;i<3;i++ )
  {
    idxsta[i] = my_sta[i];
    idxend[i] = my_end[i];
  }
  if( Allreduce( my_sta, idxsta, 3, MPI_MIN, procGrpNo ) != CPM_SUCCESS ) return false;
  if( Allreduce( my_end, idxend, 3, MPI_MAX, procGrpNo ) != CPM_SUCCESS ) return false;

  // 存在チェック
  bool bCheck = true;
  for( int i=0;i<3;i++ )
  {
    if( idxsta[i] > idxend[i] ) bCheck = false;
  }

  // 存在したとき、スタートインデクスと長さをセット
  if( bCheck )
  {
    ista = idxsta[0];
    jsta = idxsta[1];
    ksta = idxsta[2];
    ilen = idxend[0] - idxsta[0] + 1;
    jlen = idxend[1] - idxsta[1] + 1;
    klen = idxend[2] - idxsta[2] + 1;
  }

  return bCheck;
}

////////////////////////////////////////////////////////////////////////////////
// グローバルインデクスからローカルインデクスを計算
bool
cpm_ParaManager::Global2LocalIndex( int iG, int jG, int kG, int &iL, int &jL, int &kL, int procGrpNo )
{
  iL = jL = kL = 0;

  // sizeとhead
  const int *sz = GetLocalArraySize( procGrpNo );
  const int *hd = GetArrayHeadIndex( procGrpNo );
  if( !sz ) return false;

  // ローカルインデクス計算
  iL = iG - hd[0];
  jL = jG - hd[1];
  kL = kG - hd[2];

  // 内外判定
  if( iL < 0 || iL >= sz[0] ) return false;
  if( jL < 0 || jL >= sz[1] ) return false;
  if( kL < 0 || kL >= sz[2] ) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの境界が外部境界かどうかを判定
bool
cpm_ParaManager::IsOuterBoundary( cpm_FaceFlag face, int procGrpNo )
{
  // VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return false;

  // VoxelInfo
  return pVoxelInfo->IsOuterBoundary(face);
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの境界が内部境界(隣が不活性ドメイン)かどうかを判定
bool
cpm_ParaManager::IsInnerBoundary( cpm_FaceFlag face, int procGrpNo )
{
  // VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return false;

  // VoxelInfo
  return pVoxelInfo->IsInnerBoundary(face);
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファのセット
cpm_ErrorCode
cpm_ParaManager::SetBndCommBuffer( size_t maxVC, size_t maxN, int procGrpNo )
{
  if( maxVC==0 || maxN==0 )
  {
    return CPM_ERROR_BNDCOMM;
  }

  // local voxel size
  const int *sz = GetLocalArraySize(procGrpNo);
  if( !sz ) return CPM_ERROR_BNDCOMM_VOXELSIZE;

  // buffer size
  size_t nwX = size_t(sz[1]+2*maxVC) * size_t(sz[2]+2*maxVC) * maxVC * maxN;
  size_t nwY = size_t(sz[2]+2*maxVC) * size_t(sz[0]+2*maxVC) * maxVC * maxN;
  size_t nwZ = size_t(sz[0]+2*maxVC) * size_t(sz[1]+2*maxVC) * maxVC * maxN;

  // 送受信バッファ情報のインスタンス
  S_BNDCOMM_BUFFER *bufInfo = new S_BNDCOMM_BUFFER();
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_ALLOC_BUFFER;
  }

  bufInfo->m_maxVC = maxVC;
  bufInfo->m_maxN  = maxN;
  bufInfo->m_nwX   = nwX;
  bufInfo->m_nwY   = nwY;
  bufInfo->m_nwZ   = nwZ;

  // buffer
  for( int i=0;i<4;i++ )
  {
    bufInfo->m_bufX[i] = new REAL_BUF_TYPE[nwX];
    if( !bufInfo->m_bufX[i] )
    {
      delete bufInfo;
      return CPM_ERROR_BNDCOMM_ALLOC_BUFFER;
    }
    bufInfo->m_bufY[i] = new REAL_BUF_TYPE[nwY];
    if( !bufInfo->m_bufY[i] )
    {
      delete bufInfo;
      return CPM_ERROR_BNDCOMM_ALLOC_BUFFER;
    }
    bufInfo->m_bufZ[i] = new REAL_BUF_TYPE[nwZ];
    if( !bufInfo->m_bufZ[i] )
    {
      delete bufInfo;
      return CPM_ERROR_BNDCOMM_ALLOC_BUFFER;
    }
  }

  //マップにセット
  BndCommInfoMap::iterator it = m_bndCommInfoMap.find(procGrpNo);
  if( it != m_bndCommInfoMap.end() )
  {
    delete it->second;
    m_bndCommInfoMap.erase( it );
  }
  m_bndCommInfoMap.insert( std::make_pair(procGrpNo, bufInfo) );

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファサイズの取得
size_t
cpm_ParaManager::GetBndCommBufferSize( int procGrpNo )
{
  size_t mem = 0;

  // マップを検索
  if( procGrpNo < 0 )
  {
    //全体
    BndCommInfoMap::iterator it = m_bndCommInfoMap.begin();
    for( ; it!=m_bndCommInfoMap.end(); it++ )
    {
      if( !it->second ) continue;
      mem += it->second->CalcBufferSize();
    }
  }
  else
  {
    // 指定プロセスグループを検索
    S_BNDCOMM_BUFFER *bInfo = GetBndCommBuffer( procGrpNo );
    if( !bInfo ) return 0;
    mem = bInfo->CalcBufferSize();
  }

  return mem;
}
