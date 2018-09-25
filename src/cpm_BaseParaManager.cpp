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
 * @file   cpm_BaseParaManager.cpp
 * パラレルマネージャ基底クラスのソースファイル
 * @date   2012/05/31
 */
#include "cpm_BaseParaManager.h"

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_BaseParaManager::cpm_BaseParaManager()
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

  // 定義点管理マップのクリア
  m_defPointMap.clear();

}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_BaseParaManager::~cpm_BaseParaManager()
{
  // プロセスグループリストのクリア
  m_procGrpList.clear();

  // 定義点管理マップのクリア
  m_defPointMap.clear();

  // MPIのFinalize
  int flag1, flag2;
  MPI_Initialized(&flag1);
  MPI_Finalized(&flag2);
  if( flag1==true && flag2==false ) MPI_Finalize();
}

////////////////////////////////////////////////////////////////////////////////
// 初期化処理(MPI_Initは実行済みである必要がある)
cpm_ErrorCode
cpm_BaseParaManager::Initialize()
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
cpm_BaseParaManager::Initialize( int &argc, char**& argv )
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
cpm_BaseParaManager::IsParallel()
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
cpm_BaseParaManager::IsParallel() const
{
  if( m_nRank <= 1 )
  {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// プロセスグループの作成
int
cpm_BaseParaManager::CreateProcessGroup( int nproc, int *proclist, int parentProcGrpNo )
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
cpm_BaseParaManager::GetDomainType()
{
  return m_domainType;
}

////////////////////////////////////////////////////////////////////////////////
// 定義点タイプを取得
cpm_DefPointType
cpm_BaseParaManager::GetDefPointType( int procGrpNo )
{
  DefPointMap::iterator it = m_defPointMap.find(procGrpNo);
  if( it == m_defPointMap.end() ) return CPM_DEFPOINTTYPE_UNKNOWN;
  return it->second;
}

////////////////////////////////////////////////////////////////////////////////
// 全体ボクセル数を取得
const int*
cpm_BaseParaManager::GetGlobalVoxelSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalVoxelSize();
}

////////////////////////////////////////////////////////////////////////////////
// 全体頂点数を取得
const int*
cpm_BaseParaManager::GetGlobalNodeSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalNodeSize();
}

////////////////////////////////////////////////////////////////////////////////
// 全体ボクセル数または頂点数を取得
// FVMのときはボクセル数、FDMのときは頂点数を取得
const int*
cpm_BaseParaManager::GetGlobalArraySize( int procGrpNo )
{
  //定義点がVOXELのときボクセル数を取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FVM )
  {
    return GetGlobalVoxelSize(procGrpNo);
  }

  //定義点がNODEのとき頂点数を取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FDM )
  {
    return GetGlobalNodeSize(procGrpNo);
  }
  return NULL;
}

////////////////////////////////////////////////////////////////////////////////
//// 全体空間の原点を取得
const double*
cpm_BaseParaManager::GetGlobalOrigin( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 全体空間サイズを取得
const double*
cpm_BaseParaManager::GetGlobalRegion( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクのボクセル数を取得
const int*
cpm_BaseParaManager::GetLocalVoxelSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalVoxelSize();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの頂点数を取得
const int*
cpm_BaseParaManager::GetLocalNodeSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalNodeSize();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクのボクセル数または頂点数を取得
// FVMのときはボクセル数、FDMのときは頂点数を取得
const int*
cpm_BaseParaManager::GetLocalArraySize( int procGrpNo )
{
  //定義点がVOXELのときボクセル数を取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FVM ) {
    return GetLocalVoxelSize(procGrpNo);
  }

  //定義点がNODEのとき頂点数を取得
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FDM ) {
    return GetLocalNodeSize(procGrpNo);
  }

  return NULL;
}

#if 0
////////////////////////////////////////////////////////////////////////////////
// 指定面におけるレベル差を取得
int
cpm_BaseParaManager::GetNeighborLevelDiff( cpm_FaceFlag face, int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return 0;

  return pVoxelInfo->GetNeighborLevelDiff( face );
}
#endif

////////////////////////////////////////////////////////////////////////////////
// パディングサイズ取得処理(静的関数)
// 暫定対応として、偶数のときに+1するものとする
int
cpm_BaseParaManager::GetPaddingSize1D( const int size, const int vc )
{
  int npad = 0;
  if( (size+2*vc)%2 == 0 )
  {
    npad = 1;
  }
  return npad;
}

#define _ALL_DIM_PAD_
////////////////////////////////////////////////////////////////////////////////
// パディングサイズ取得処理(静的関数)
// 暫定対応として、1次元目(最内ループ)のみを偶数のときに+1するものとする
void
cpm_BaseParaManager::GetPaddingSize( CPM_ARRAY_SHAPE atype, const int *size, const int vc, int *pad_size, int nmax )
{
  if( atype==CPM_ARRAY_V3D || atype==CPM_ARRAY_V3DEX )
  {
    nmax = 3;
  }

  switch( atype )
  {
  case CPM_ARRAY_S3D:
    pad_size[0] = GetPaddingSize1D( size[0], vc );
    pad_size[1] = 0;
    pad_size[2] = 0;
#ifdef _ALL_DIM_PAD_
    pad_size[1] = GetPaddingSize1D( size[1], vc );
    pad_size[2] = GetPaddingSize1D( size[2], vc );
#endif
    break;
  case CPM_ARRAY_V3D:
  case CPM_ARRAY_S4D:
    pad_size[0] = GetPaddingSize1D( size[0], vc );
    pad_size[1] = 0;
    pad_size[2] = 0;
    pad_size[3] = 0;
#ifdef _ALL_DIM_PAD_
    pad_size[1] = GetPaddingSize1D( size[1], vc );
    pad_size[2] = GetPaddingSize1D( size[2], vc );
    pad_size[3] = GetPaddingSize1D( nmax, 0 );
#endif
    break;
  case CPM_ARRAY_V3DEX:
  case CPM_ARRAY_S4DEX:
    pad_size[0] = GetPaddingSize1D( nmax, 0 );
    pad_size[1] = 0;
    pad_size[2] = 0;
    pad_size[3] = 0;
#ifdef _ALL_DIM_PAD_
    pad_size[1] = GetPaddingSize1D( size[0], vc );
    pad_size[2] = GetPaddingSize1D( size[1], vc );
    pad_size[3] = GetPaddingSize1D( size[2], vc );
#endif
    break;
  }
}

#ifdef _DEBUG
#include <fstream>
// debug write
void
cpm_BaseParaManager::printVoxelInfo(int myrank)
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
          const int    *gvox = pV->GetGlobalVoxelSize();
          const int    *gnod = pV->GetGlobalNodeSize();
          cpm_DefPointType dtype = GetDefPointType();
          const int    *gary = pV->GetGlobalArraySize(dtype);
          std::cout << "  +----------------------------------" << std::endl;
          std::cout << "  global div  = " << gdiv[0] << "," << gdiv[1] << "," << gdiv[2] << std::endl;
          std::cout << "  global org  = " << gorg[0] << "," << gorg[1] << "," << gorg[2] << std::endl;
          std::cout << "  global pch  = " << gpch[0] << "," << gpch[1] << "," << gpch[2] << std::endl;
          std::cout << "  global rgn  = " << grgn[0] << "," << grgn[1] << "," << grgn[2] << std::endl;
          std::cout << "  global vox  = " << gvox[0] << "," << gvox[1] << "," << gvox[2] << std::endl;
          std::cout << "  global nod  = " << gnod[0] << "," << gnod[1] << "," << gnod[2] << std::endl;
          std::cout << "  global array= " << gary[0] << "," << gary[1] << "," << gary[2] << std::endl;
          std::cout << "  DefPointType= " << dtype   << std::endl;

          // ローカル空間
          const double *lorg = pV->GetLocalOrigin();
          const double *lpch = pV->GetPitch();
          const double *lrgn = pV->GetLocalRegion();
          const int    *lvox = pV->GetLocalVoxelSize();
          const int    *lnod = pV->GetLocalNodeSize();
          const int    *lary = pV->GetLocalArraySize(dtype);
          const int    *lpos = pV->GetDivPos();
          const int    *head = pV->GetVoxelHeadIndex();
          const int    *tail = pV->GetVoxelTailIndex();
          const int    *nhead = pV->GetNodeHeadIndex();
          const int    *ntail = pV->GetNodeTailIndex();
          const int    *ahead = pV->GetArrayHeadIndex(dtype);
          const int    *atail = pV->GetArrayTailIndex(dtype);
//          const int    *neig = pV->GetNeighborRankID();
//          const int    *peri = pV->GetPeriodicRankID();
          std::cout << "  +----------------------------------" << std::endl;
          std::cout << "  local  org  = " << lorg[0] << "," << lorg[1] << "," << lorg[2] << std::endl;
          std::cout << "  local  pch  = " << lpch[0] << "," << lpch[1] << "," << lpch[2] << std::endl;
          std::cout << "  local  rgn  = " << lrgn[0] << "," << lrgn[1] << "," << lrgn[2] << std::endl;
          std::cout << "  local  vox  = " << lvox[0] << "," << lvox[1] << "," << lvox[2] << std::endl;
          std::cout << "  local  nod  = " << lnod[0] << "," << lnod[1] << "," << lnod[2] << std::endl;
          std::cout << "  local  array= " << lary[0] << "," << lary[1] << "," << lary[2] << std::endl;
          std::cout << "  local  pos  = " << lpos[0] << "," << lpos[1] << "," << lpos[2] << std::endl;
//        std::cout << "  local  head = " << head[0] << "," << head[1] << "," << head[2] << std::endl;
//        std::cout << "  local  tail = " << tail[0] << "," << tail[1] << "," << tail[2] << std::endl;
          std::cout << "  local  vox head  = " <<  head[0] << "," <<  head[1] << "," <<  head[2] << std::endl;
          std::cout << "  local  vox tail  = " <<  tail[0] << "," <<  tail[1] << "," <<  tail[2] << std::endl;
          std::cout << "  local  nod head  = " << nhead[0] << "," << nhead[1] << "," << nhead[2] << std::endl;
          std::cout << "  local  nod tail  = " << ntail[0] << "," << ntail[1] << "," << ntail[2] << std::endl;
          std::cout << "  local  array head= " << ahead[0] << "," << ahead[1] << "," << ahead[2] << std::endl;
          std::cout << "  local  array tail= " << atail[0] << "," << atail[1] << "," << atail[2] << std::endl;
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
      const int    *gvox = pV->GetGlobalVoxelSize();
      const int    *gnod = pV->GetGlobalNodeSize();
      cpm_DefPointType dtype = GetDefPointType();
      const int    *gary = pV->GetGlobalArraySize(dtype);
      ofs << "  +----------------------------------" << std::endl;
      ofs << "  global div  = " << gdiv[0] << "," << gdiv[1] << "," << gdiv[2] << std::endl;
      ofs << "  global org  = " << gorg[0] << "," << gorg[1] << "," << gorg[2] << std::endl;
      ofs << "  global pch  = " << gpch[0] << "," << gpch[1] << "," << gpch[2] << std::endl;
      ofs << "  global rgn  = " << grgn[0] << "," << grgn[1] << "," << grgn[2] << std::endl;
      ofs << "  global vox  = " << gvox[0] << "," << gvox[1] << "," << gvox[2] << std::endl;
      ofs << "  global nod  = " << gnod[0] << "," << gnod[1] << "," << gnod[2] << std::endl;
      ofs << "  global array= " << gary[0] << "," << gary[1] << "," << gary[2] << std::endl;
      ofs << "  DefPointType= " << dtype   << std::endl;

      // ローカル空間
      const double *lorg = pV->GetLocalOrigin();
      const double *lpch = pV->GetPitch();
      const double *lrgn = pV->GetLocalRegion();
      const int    *lvox = pV->GetLocalVoxelSize();
      const int    *lnod = pV->GetLocalNodeSize();
      const int    *lary = pV->GetLocalArraySize(dtype);
      const int    *lpos = pV->GetDivPos();
      const int    *head = pV->GetVoxelHeadIndex();
      const int    *tail = pV->GetVoxelTailIndex();
      const int    *nhead = pV->GetNodeHeadIndex();
      const int    *ntail = pV->GetNodeTailIndex();
      const int    *ahead = pV->GetArrayHeadIndex(dtype);
      const int    *atail = pV->GetArrayTailIndex(dtype);
//      const int    *neig = pV->GetNeighborRankID();
//      const int    *peri = pV->GetPeriodicRankID();
      ofs << "  +----------------------------------" << std::endl;
      ofs << "  local  org  = " << lorg[0] << "," << lorg[1] << "," << lorg[2] << std::endl;
      ofs << "  local  pch  = " << lpch[0] << "," << lpch[1] << "," << lpch[2] << std::endl;
      ofs << "  local  rgn  = " << lrgn[0] << "," << lrgn[1] << "," << lrgn[2] << std::endl;
      ofs << "  local  vox  = " << lvox[0] << "," << lvox[1] << "," << lvox[2] << std::endl;
      ofs << "  local  nod  = " << lnod[0] << "," << lnod[1] << "," << lnod[2] << std::endl;
      ofs << "  local  array= " << lary[0] << "," << lary[1] << "," << lary[2] << std::endl;
      ofs << "  local  pos = " << lpos[0] << "," << lpos[1] << "," << lpos[2] << std::endl;
//    ofs << "  local  head= " << head[0] << "," << head[1] << "," << head[2] << std::endl;
//    ofs << "  local  tail= " << tail[0] << "," << tail[1] << "," << tail[2] << std::endl;
      ofs << "  local  vox head=  " <<  head[0] << "," <<  head[1] << "," <<  head[2] << std::endl;
      ofs << "  local  vox tail=  " <<  tail[0] << "," <<  tail[1] << "," <<  tail[2] << std::endl;
      ofs << "  local  nod head=  " << nhead[0] << "," << nhead[1] << "," << nhead[2] << std::endl;
      ofs << "  local  nod tail=  " << ntail[0] << "," << ntail[1] << "," << ntail[2] << std::endl;
      ofs << "  local  array head=" << ahead[0] << "," << ahead[1] << "," << ahead[2] << std::endl;
      ofs << "  local  array tail=" << atail[0] << "," << atail[1] << "," << atail[2] << std::endl;
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
