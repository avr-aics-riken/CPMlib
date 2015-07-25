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
 * @file   cpm_ParaManager.cpp
 * パラレルマネージャクラスのソースファイル
 * @date   2012/05/31
 */
#include "cpm_ParaManager.h"
#include "cpm_ParaManagerCART.h"
#include "cpm_ParaManagerLMR.h"

/** 並列管理クラスの自動破棄管理クラス
 */
class C_PARAMANAGER
{
friend class cpm_ParaManager;
private:
  cpm_ParaManager *pParaManager;
  C_PARAMANAGER()
  {
    pParaManager = NULL;
  }
  ~C_PARAMANAGER()
  {
    delete pParaManager;
    pParaManager = NULL;
  }
  static C_PARAMANAGER* get_instance()
  {
    static C_PARAMANAGER instance;
    return &instance;
  }
};

////////////////////////////////////////////////////////////////////////////////
// 唯一のインスタンスの取得
cpm_ParaManager*
cpm_ParaManager::get_instance()
{
  // 管理クラスのインスタンスを取得
  C_PARAMANAGER *instance = C_PARAMANAGER::get_instance();
  if( !instance ) return NULL;

  // ポインタ
  return instance->pParaManager;
}

////////////////////////////////////////////////////////////////////////////////
// 唯一のインスタンスの取得(initialize処理も実行)
cpm_ParaManager*
cpm_ParaManager::get_instance(int &argc, char**& argv, cpm_DomainType domainType)
{
  // 管理クラスのインスタンスを取得
  C_PARAMANAGER *instance = C_PARAMANAGER::get_instance();
  if( !instance ) return NULL;

  // インスタンス済みの場合、ポインタを返す
  if( instance->pParaManager )
  {
    return instance->pParaManager;
  }

  // インスタンス
  if( domainType == CPM_DOMAIN_CARTESIAN )
  {
    instance->pParaManager = new cpm_ParaManagerCART();
  }
  else if( domainType == CPM_DOMAIN_LMR )
  {
    instance->pParaManager = new cpm_ParaManagerLMR();
  }
  else
  {
    return NULL;
  }

  // MPI_Init
  if( instance->pParaManager->Initialize(argc,argv) != CPM_SUCCESS )
  {
    return NULL;
  }

  // ポインタ
  return instance->pParaManager;
}

////////////////////////////////////////////////////////////////////////////////
// 唯一のインスタンスの取得(initialize処理も実行)
// fortranインターフェイス用
cpm_ParaManager*
cpm_ParaManager::get_instance(cpm_DomainType domainType)
{
  // 管理クラスのインスタンスを取得
  C_PARAMANAGER *instance = C_PARAMANAGER::get_instance();
  if( !instance ) return NULL;

  // インスタンス済みの場合、ポインタを返す
  if( instance->pParaManager )
  {
    return instance->pParaManager;
  }

  // インスタンス
  if( domainType == CPM_DOMAIN_CARTESIAN )
  {
    instance->pParaManager = new cpm_ParaManagerCART();
  }
  else if( domainType == CPM_DOMAIN_LMR )
  {
    instance->pParaManager = new cpm_ParaManagerLMR();
  }
  else
  {
    return NULL;
  }

  // MPI_Init
  if( instance->pParaManager->Initialize() != CPM_SUCCESS )
  {
    return NULL;
  }

  // ポインタ
  return instance->pParaManager;
}

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_ParaManager::cpm_ParaManager()
  : cpm_Base()
{
  // 並列数、ランク番号
  m_nRank  = 1;
  m_rankNo = 0;

  // 領域分割タイプ
  m_domainType = CPM_DOMAIN_UNKNOWN;

  // プロセスグループリストをクリア
  m_procGrpList.clear();
  m_procGrpList.push_back(MPI_COMM_WORLD);

  // 領域管理情報マップのクリア
  m_voxelInfoMap.clear();
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_ParaManager::~cpm_ParaManager()
{
  // プロセスグループリストのクリア
  m_procGrpList.clear();

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

  // MPIのFinalize
  int flag1, flag2;
  MPI_Initialized(&flag1);
  MPI_Finalized(&flag2);
  if( flag1==true && flag2==false ) MPI_Finalize();
}

////////////////////////////////////////////////////////////////////////////////
// 初期化処理(MPI_Initは実行済みである必要がある)
cpm_ErrorCode
cpm_ParaManager::Initialize()
{
  m_nRank  = 1;
  m_rankNo = 0;

  // MPI_Init実行済みかチェック
  int flag1;
  MPI_Initialized(&flag1);
  if( flag1==false )
  {
    return CPM_ERROR_NO_MPI_INIT;
  }

  // プロセス並列数の取得
  if( MPI_Comm_size(MPI_COMM_WORLD, &m_nRank) != MPI_SUCCESS )
  {
    std::cerr << "MPI_Comm_size error." << std::endl;
    return CPM_ERROR_MPI;
  }
  if( m_nRank < 1 ) m_nRank = 1;

  // 自ランク番号の取得
  if( MPI_Comm_rank(MPI_COMM_WORLD, &m_rankNo) != MPI_SUCCESS )
  {
    std::cerr << "MPI_Comm_rank error." << std::endl;
    return CPM_ERROR_MPI;
  }
  if( m_rankNo < 0 ) m_rankNo = 0;

#ifdef _DEBUG
  Barrier();
  std::cout << std::flush;
  Barrier();
  for( int i=0;i<m_nRank;i++ )
  {
    if( i==m_rankNo )
      std::cout << "[" << m_rankNo << "] #rank=" << m_nRank << " parallel=" << IsParallel() << std::endl;
    std::cout << std::flush;
    Barrier();
  }
  std::cout << std::flush;
  Barrier();
#endif

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 初期化処理(MPI_Initも実行)
cpm_ErrorCode
cpm_ParaManager::Initialize( int &argc, char**& argv )
{
  m_nRank  = 1;
  m_rankNo = 0;

  // MPI_Init
  int flag1;
  MPI_Initialized(&flag1);
  if( flag1==false )
  {
    if( MPI_Init(&argc,&argv) != MPI_SUCCESS )
    {
      std::cerr << "MPI_Init error." << std::endl;
      return CPM_ERROR_MPI;
    }
  }

  // initialize
  return Initialize();
}

////////////////////////////////////////////////////////////////////////////////
// 並列実行であるかチェックする
bool
cpm_ParaManager::IsParallel()
{
  if( m_nRank <= 1 )
  {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 並列実行であるかチェックする(const)
bool
cpm_ParaManager::IsParallel() const
{
  if( m_nRank <= 1 )
  {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// カーテシアン用の領域分割
cpm_ErrorCode
cpm_ParaManager::VoxelInit( cpm_GlobalDomainInfo* domainInfo, size_t maxVC, size_t maxN, int procGrpNo )
{
  return CPM_ERROR_DOMAINTYPE_VOXELINIT;
}

////////////////////////////////////////////////////////////////////////////////
// カーテシアン用の領域分割(領域分割数を指定)
cpm_ErrorCode
cpm_ParaManager::VoxelInit( int div[3], int vox[3], double origin[3], double region[3]
                          , size_t maxVC, size_t maxN, cpm_DivPolicy divPolicy, int procGrpNo )
{
  return CPM_ERROR_DOMAINTYPE_VOXELINIT;
}

////////////////////////////////////////////////////////////////////////////////
// カーテシアン用の領域分割(プロセスグループのランク数で自動領域分割)
cpm_ErrorCode
cpm_ParaManager::VoxelInit( int vox[3], double origin[3], double pitch[3]
                          , size_t maxVC, size_t maxN, cpm_DivPolicy divPolicy
                          , int procGrpNo )
{
  return CPM_ERROR_DOMAINTYPE_VOXELINIT;
}

////////////////////////////////////////////////////////////////////////////////
// カーテシアン用の領域分割(ActiveSubdomainファイル、領域分割数を指定)
cpm_ErrorCode
cpm_ParaManager::VoxelInit_Subdomain( int div[3], int vox[3], double origin[3], double region[3]
                                    , std::string subDomainFile
                                    , size_t maxVC, size_t maxN, int procGrpNo )
{
  return CPM_ERROR_DOMAINTYPE_VOXELINIT;
}

////////////////////////////////////////////////////////////////////////////////
// カーテシアン用の領域分割(ActiveSubdomainファイルを指定)
cpm_ErrorCode
cpm_ParaManager::VoxelInit_Subdomain( int vox[3], double origin[3], double pitch[3]
                                    , std::string subDomainFile
                                    , size_t maxVC, size_t maxN, int procGrpNo )
{
  return CPM_ERROR_DOMAINTYPE_VOXELINIT;
}

////////////////////////////////////////////////////////////////////////////////
// LMR用の領域分割
cpm_ErrorCode
cpm_ParaManager::VoxelInit_LMR( std::string treeFile
                              , size_t maxVC, size_t maxN, int procGrpNo )
{
  return CPM_ERROR_DOMAINTYPE_VOXELINIT;
}

////////////////////////////////////////////////////////////////////////////////
// プロセスグループの作成
int
cpm_ParaManager::CreateProcessGroup( int nproc, int *proclist, int parentProcGrpNo )
{
  // 入力チェック
  if( parentProcGrpNo < 0 || nproc < 1 || !proclist )
  {
    return -1;
  }

  // 親のMPIコミュニケータを取得
  MPI_Comm parentComm = GetMPI_Comm( parentProcGrpNo );
  if( IsCommNull( parentComm ) )
  {
    return -1;
  }

  // 親のグループを取得
  MPI_Group parentGrp;
  MPI_Comm_group(parentComm, &parentGrp);

  // 新しいグループを生成
  MPI_Group newGrp;
  MPI_Group_incl(parentGrp, nproc, proclist, &newGrp);

  // 新しいコミュニケータを生成
  MPI_Comm  newComm;
  MPI_Comm_create(parentComm, newGrp, &newComm);
  if( IsCommNull(newComm) )
  {
    // 時ランクが含まれない
    return -1;
  }

  // コミュニケータを登録
  m_procGrpList.push_back(newComm);

  // プロセスグループの番号
  return m_procGrpList.size()-1;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割タイプを取得
cpm_DomainType
cpm_ParaManager::GetDomainType()
{
  return m_domainType;
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
// 全体ボクセル数を取得
const int*
cpm_ParaManager::GetGlobalVoxelSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalVoxelSize();
}

////////////////////////////////////////////////////////////////////////////////
// 全体空間の原点を取得
const double*
cpm_ParaManager::GetGlobalOrigin( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 全体空間サイズを取得
const double*
cpm_ParaManager::GetGlobalRegion( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクのボクセル数を取得
const int*
cpm_ParaManager::GetLocalVoxelSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalVoxelSize();
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
// 指定面における自ランクの隣接ランク番号を取得
const int*
cpm_ParaManager::GetNeighborRankList( cpm_FaceFlag face, int &num, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNeighborRankList( face, num );
}

////////////////////////////////////////////////////////////////////////////////
// 指定面における自ランクの周期境界の隣接ランク番号を取得
const int*
cpm_ParaManager::GetPeriodicRankList( cpm_FaceFlag face, int &num, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPeriodicRankList( face, num );
}

////////////////////////////////////////////////////////////////////////////////
// 指定面におけるレベル差を取得
int
cpm_ParaManager::GetNeighborLevelDiff( cpm_FaceFlag face, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return 0;

  return pVoxelInfo->GetNeighborLevelDiff( face );
}

////////////////////////////////////////////////////////////////////////////////
// 指定idを含む全体ボクセル空間のインデクス範囲を取得
bool
cpm_ParaManager::GetBndIndexExtGc( int id, int *array, int vc
                                 , int &ista, int &jsta, int &ksta
                                 , int &ilen, int &jlen, int &klen
                                 , int procGrpNo )
{
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// 指定idを含む全体ボクセル空間のインデクス範囲を取得
bool
cpm_ParaManager::GetBndIndexExtGc( int id, int *array
                               , int imax, int jmax, int kmax, int vc
                               , int &ista, int &jsta, int &ksta
                               , int &ilen, int &jlen, int &klen
                               , int procGrpNo )
{
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの境界が外部境界かどうかを判定
bool
cpm_ParaManager::Global2LocalIndex( int iG, int jG, int kG, int &iL, int &jL, int &kL, int procGrpNo )
{
  iL = jL = kL = 0;

  // sizeとhead
  const int *sz = GetLocalVoxelSize( procGrpNo );
  const int *hd = GetVoxelHeadIndex( procGrpNo );
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
  return CPM_ERROR_DOMAINTYPE_SETBNDCOMMBUF;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファサイズの取得
size_t
cpm_ParaManager::GetBndCommBufferSize( int procGrpNo )
{
  return size_t(0);
}

#ifdef _DEBUG
#include <fstream>
// debug write
void
cpm_ParaManager::printVoxelInfo(int myrank)
{
  if( myrank < 0 )
  {
    Barrier();
    std::cout << std::flush;
    Barrier();
    for( int i=0;i<m_nRank;i++ )
    {
      if( i==m_rankNo )
      {
        std::cout << "####### VoxelInfo [" << i << " @ MPI_COMM_WORLD] #######" << std::endl;
        VoxelInfoMap::iterator itv = m_voxelInfoMap.begin();
        for( ;itv!=m_voxelInfoMap.end();itv++ )
        {
          // プロセスグループ番号、ボクセル情報、ランク番号
          int procGrpNo = itv->first;
          cpm_VoxelInfo *pV = itv->second;
          if( !pV ) continue;
          int rankNo = GetMyRankID(procGrpNo);
          std::cout << " *process group No [" << procGrpNo << "], rankNo[" << rankNo << "]" << std::endl;

          // 全体空間
          const int    *gdiv = pV->GetDivNum();
          const double *gorg = pV->GetGlobalOrigin();
          const double *gpch = pV->GetGlobalPitch();
          const double *grgn = pV->GetGlobalRegion();
          const int       *gvox = pV->GetGlobalVoxelSize();
          std::cout << "  +----------------------------------" << std::endl;
          std::cout << "  global div = " << gdiv[0] << "," << gdiv[1] << "," << gdiv[2] << std::endl;
          std::cout << "  global org = " << gorg[0] << "," << gorg[1] << "," << gorg[2] << std::endl;
          std::cout << "  global pch = " << gpch[0] << "," << gpch[1] << "," << gpch[2] << std::endl;
          std::cout << "  global rgn = " << grgn[0] << "," << grgn[1] << "," << grgn[2] << std::endl;
          std::cout << "  global vox = " << gvox[0] << "," << gvox[1] << "," << gvox[2] << std::endl;

          // ローカル空間
          const double *lorg = pV->GetLocalOrigin();
          const double *lpch = pV->GetPitch();
          const double *lrgn = pV->GetLocalRegion();
          const int    *lvox = pV->GetLocalVoxelSize();
          const int    *lpos = pV->GetDivPos();
          const int    *head = pV->GetVoxelHeadIndex();
          const int    *tail = pV->GetVoxelTailIndex();
//          const int    *neig = pV->GetNeighborRankID();
//          const int    *peri = pV->GetPeriodicRankID();
          std::cout << "  +----------------------------------" << std::endl;
          std::cout << "  local  org = " << lorg[0] << "," << lorg[1] << "," << lorg[2] << std::endl;
          std::cout << "  local  pch = " << lpch[0] << "," << lpch[1] << "," << lpch[2] << std::endl;
          std::cout << "  local  rgn = " << lrgn[0] << "," << lrgn[1] << "," << lrgn[2] << std::endl;
          std::cout << "  local  vox = " << lvox[0] << "," << lvox[1] << "," << lvox[2] << std::endl;
          std::cout << "  local  pos = " << lpos[0] << "," << lpos[1] << "," << lpos[2] << std::endl;
          std::cout << "  local  head= " << head[0] << "," << head[1] << "," << head[2] << std::endl;
          std::cout << "  local  tail= " << tail[0] << "," << tail[1] << "," << tail[2] << std::endl;
//          std::cout << "  local  neig= " << neig[0] << "," << neig[1] << "," << neig[2] << ","
//                                         << neig[3] << "," << neig[4] << "," << neig[5] << std::endl;
//          std::cout << "  local  peri= " << peri[0] << "," << peri[1] << "," << peri[2] << ","
//                                         << peri[3] << "," << peri[4] << "," << peri[5] << std::endl;
          const char fname[6][3] = {"-X", "+X", "-Y", "+Y", "-Z", "+Z"};
          for( int f=0;f<6;f++ )
          {
            int diff = pV->GetNeighborLevelDiff( cpm_FaceFlag(f) );
            std::cout << "  local  diff " << fname[f] << "=" << diff << std::endl;
          }
          for( int f=0;f<6;f++ )
          {
            int num = 0;
            const int *neig = pV->GetNeighborRankList( cpm_FaceFlag(f), num );
            std::cout << "  local  neig " << fname[f] << "=";
            for( int m=0;m<num;m++ )
            {
              std::cout << neig[m] << ",";
            }
            std::cout << std::endl;
          }
          for( int f=0;f<6;f++ )
          {
            int num = 0;
            const int *peri = pV->GetPeriodicRankList( cpm_FaceFlag(f), num );
            std::cout << "  local  peri " << fname[f] << "=";
            for( int m=0;m<num;m++ )
            {
              std::cout << peri[m] << ",";
            }
            std::cout << std::endl;
          }
        }
      }
      Barrier();
      std::cout << std::flush;
      Barrier();
    }
    Barrier();
    std::cout << std::flush;
    Barrier();
  }
  else
  {
    char fname[512];
    sprintf ( fname, "vinfo_%04d.log", myrank );
    std::ofstream ofs( fname );
    ofs << "####### VoxelInfo #######" << std::endl;
    VoxelInfoMap::iterator itv = m_voxelInfoMap.begin();
    for( ;itv!=m_voxelInfoMap.end();itv++ )
    {
      // プロセスグループ番号、ボクセル情報、ランク番号
      int procGrpNo = itv->first;
      cpm_VoxelInfo *pV = itv->second;
      if( !pV ) continue;
      int rankNo = GetMyRankID(procGrpNo);
      ofs << " *process group No [" << procGrpNo << "], rankNo[" << rankNo << "]" << std::endl;

      // 全体空間
      const int    *gdiv = pV->GetDivNum();
      const double *gorg = pV->GetGlobalOrigin();
      const double *gpch = pV->GetGlobalPitch();
      const double *grgn = pV->GetGlobalRegion();
      const int       *gvox = pV->GetGlobalVoxelSize();
      ofs << "  +----------------------------------" << std::endl;
      ofs << "  global div = " << gdiv[0] << "," << gdiv[1] << "," << gdiv[2] << std::endl;
      ofs << "  global org = " << gorg[0] << "," << gorg[1] << "," << gorg[2] << std::endl;
      ofs << "  global pch = " << gpch[0] << "," << gpch[1] << "," << gpch[2] << std::endl;
      ofs << "  global rgn = " << grgn[0] << "," << grgn[1] << "," << grgn[2] << std::endl;
      ofs << "  global vox = " << gvox[0] << "," << gvox[1] << "," << gvox[2] << std::endl;

      // ローカル空間
      const double *lorg = pV->GetLocalOrigin();
      const double *lpch = pV->GetPitch();
      const double *lrgn = pV->GetLocalRegion();
      const int    *lvox = pV->GetLocalVoxelSize();
      const int    *lpos = pV->GetDivPos();
      const int    *head = pV->GetVoxelHeadIndex();
      const int    *tail = pV->GetVoxelTailIndex();
//      const int    *neig = pV->GetNeighborRankID();
//      const int    *peri = pV->GetPeriodicRankID();
      ofs << "  +----------------------------------" << std::endl;
      ofs << "  local  org = " << lorg[0] << "," << lorg[1] << "," << lorg[2] << std::endl;
      ofs << "  local  pch = " << lpch[0] << "," << lpch[1] << "," << lpch[2] << std::endl;
      ofs << "  local  rgn = " << lrgn[0] << "," << lrgn[1] << "," << lrgn[2] << std::endl;
      ofs << "  local  vox = " << lvox[0] << "," << lvox[1] << "," << lvox[2] << std::endl;
      ofs << "  local  pos = " << lpos[0] << "," << lpos[1] << "," << lpos[2] << std::endl;
      ofs << "  local  head= " << head[0] << "," << head[1] << "," << head[2] << std::endl;
      ofs << "  local  tail= " << tail[0] << "," << tail[1] << "," << tail[2] << std::endl;
//      ofs << "  local  neig= " << neig[0] << "," << neig[1] << "," << neig[2] << ","
//                               << neig[3] << "," << neig[4] << "," << neig[5] << std::endl;
//      ofs << "  local  peri= " << peri[0] << "," << peri[1] << "," << peri[2] << ","
//                               << peri[3] << "," << peri[4] << "," << peri[5] << std::endl;
      const char fname[6][3] = {"-X", "+X", "-Y", "+Y", "-Z", "+Z"};
      for( int f=0;f<6;f++ )
      {
        int diff = pV->GetNeighborLevelDiff( cpm_FaceFlag(f) );
        ofs << "  local  diff " << fname[f] << "=" << diff << std::endl;
      }
      for( int f=0;f<6;f++ )
      {
        int num = 0;
        const int *neig = pV->GetNeighborRankList( cpm_FaceFlag(f), num );
        ofs << "  local  neig " << fname[f] << "=";
        for( int m=0;m<num;m++ )
        {
          ofs << neig[m] << ",";
        }
        ofs << std::endl;
      }
      for( int f=0;f<6;f++ )
      {
        int num = 0;
        const int *peri = pV->GetPeriodicRankList( cpm_FaceFlag(f), num );
        ofs << "  local  peri " << fname[f] << "=";
        for( int m=0;m<num;m++ )
        {
          ofs << peri[m] << ",";
        }
        ofs << std::endl;
      }
    }

    ofs << "####### rank map #######" << std::endl;
    itv = m_voxelInfoMap.begin();
    for( ;itv!=m_voxelInfoMap.end();itv++ )
    {
      // プロセスグループ番号、ボクセル情報、ランク番号
      int procGrpNo = itv->first;
      cpm_VoxelInfo *pV = itv->second;
      if( !pV ) continue;
      int rankNo = GetMyRankID(procGrpNo);
      ofs << " *process group No [" << procGrpNo << "], rankNo[" << rankNo << "]" << std::endl;
#if 0
      // マップ、領域分割数
      int *map = pV->m_rankMap;
      const int *div = pV->GetDivNum();
      ofs << " *div = " << div[0] << " x " << div[1] << " x " << div[2] << std::endl;
      ofs << " *map-------------------------------" << std::endl;
      for( int k=div[2]-1;k>=0;k-- )
      {
        ofs << "  k=" << k << std::endl;
        for( int j=div[1]-1;j>=0;j-- )
        {
           ofs << "    ";
           for( int i=0;i<div[0];i++ )
           {
             ofs.width(3);
             ofs.setf(std::ios::left, std::ios::adjustfield);
             ofs << map[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)];
             if( i < div[0]-1 ) ofs << ",";
           }
           ofs << std::endl;
        }
      }
#endif
    }
    ofs.close();
  }
}
#endif
