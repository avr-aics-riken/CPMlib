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
 * @file   cpm_VoxelInfoLMR.cpp
 * VOXEL空間情報クラスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "cpm_VoxelInfoLMR.h"
using namespace BCMFileIO;

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_VoxelInfoLMR::cpm_VoxelInfoLMR()
  : cpm_VoxelInfo()
{
  m_octree = NULL;
  m_neighborInfo = NULL;
  for( int m=0;m<6;m++ )
  {
    for( int n=0;n<4;n++ )
    {
      m_neighborRankID_LMR[m][n] = getRankNull();
      m_periodicRankID_LMR[m][n] = getRankNull();
    }
    m_neighborLevelDiff[m] = 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_VoxelInfoLMR::~cpm_VoxelInfoLMR()
{
  delete m_octree;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割情報の生成
cpm_ErrorCode
cpm_VoxelInfoLMR::Init( MPI_Comm comm, std::string treeFile )
{
  cpm_ErrorCode ret = CPM_SUCCESS;

  // 入力チェック
  if( IsCommNull(comm) )
  {
    return CPM_ERROR_MPI_INVALID_COMM;
  }

  // 入力をコピー
  m_comm = comm;

  // ランク数、ランク番号をセット
  MPI_Comm_size(m_comm, &m_nRank);
  MPI_Comm_rank(m_comm, &m_rankNo);

  // 領域情報の読み込み
  S_OCT_DOMAIN_INFO domainInfo;
  if( (ret = cpm_TextParserDomainLMR::Read(treeFile, domainInfo)) != CPM_SUCCESS )
  {
    return ret;
  }

  // 木情報を読み込み
  std::vector<Pedigree> pedigrees;
  if( (ret = LoadOctreeFile(domainInfo.octFile, m_octHeader, pedigrees)) != CPM_SUCCESS )
  {
    return ret;
  }

#if 1
if( m_rankNo==0 )
{
    domainInfo.print();
    std::cout << "*** header info" << std::endl;
    std::cout << "org = " << m_octHeader.org[0] << "," << m_octHeader.org[1] << "," << m_octHeader.org[2] << std::endl;
    std::cout << "rgn = " << m_octHeader.rgn[0] << "," << m_octHeader.rgn[1] << "," << m_octHeader.rgn[2] << std::endl;
    std::cout << "rootDim = " << m_octHeader.rootDims[0] << "," << m_octHeader.rootDims[1] << "," << m_octHeader.rootDims[2] << std::endl;
    std::cout << "maxLevel = " << m_octHeader.maxLevel << std::endl;
    std::cout << "numLeaf = " << m_octHeader.numLeaf << std::endl;
    std::cout << "padding = " << m_octHeader.padding << std::endl;
    std::cout << "*** pedigree info" << std::endl;
    std::cout << "num pedigree=" << pedigrees.size() << std::endl;
    for( size_t i=0;i<pedigrees.size();i++ )
    {
      std::cout << pedigrees[i] << std::endl;
    }
  }
#endif

  // 並列数=リーフ数のチェック
  if( m_nRank != m_octHeader.numLeaf )
  {
    return CPM_ERROR_LMR_MISMATCH_NP_NUMLEAF;
  }

  // RootGrid、BCMOctreeの生成:
  RootGrid *rootGrid = new RootGrid(m_octHeader.rootDims[0], m_octHeader.rootDims[1], m_octHeader.rootDims[2]);
  m_octree  = new BCMOctree(rootGrid, pedigrees);
  if( !m_octree )
  {
    return CPM_ERROR_LMR_INVALID_OCTFILE;
  }

  // 自身の担当リーフを決定
  std::vector<Node*> &leafNodeArray = m_octree->getLeafNodeArray();
  m_node = leafNodeArray[m_rankNo];
#if 1
if( m_rankNo==0 )
{
  std::cout << "*** pedigree info @ octree" << std::endl;
  for( size_t i=0;i<leafNodeArray.size();i++ )
  {
    int BlockID = leafNodeArray[i]->getBlockID();
    Vec3d org = m_octree->getOrigin(leafNodeArray[i]);
    std::cout << leafNodeArray[i]->getPedigree() << " , " 
              << "org=(" <<org[0] << "," << org[1] << "," << org[2] << ")" << " , "
              << "blocID=" << BlockID
              << std::endl;
  }
}
std::cout << "*** Node @ " << m_rankNo << m_node->getPedigree() << std::endl;
#endif

  // 領域情報のセット
  SetGlobaliDomainInfo( domainInfo );
  SetLocalDomainInfo( domainInfo );

  // 隣接情報の取得
  SetNeighborInfo();






  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// 木情報ファイルの読み込み
cpm_ErrorCode
cpm_VoxelInfoLMR::LoadOctreeFile( std::string octFile, OctHeader &header, std::vector<Pedigree> &pedigrees )
{
  cpm_ErrorCode ret = CPM_SUCCESS;
  bool isNeedSwap = false;

  // ファイルオープン
  FILE *fp = fopen( octFile.c_str(), "rb" );
  if( !fp )
  {
    return CPM_ERROR_LMR_OPEN_OCTFILE;
  }

  // ヘッダー読み込み
  if( (ret = LoadOctreeHeader(fp, header, isNeedSwap)) != CPM_SUCCESS )
  {
    fclose(fp);
    return ret;
  }

  // ぺディグリー読み込み
  pedigrees.resize(header.numLeaf);
  fread(&pedigrees[0], sizeof(Pedigree), header.numLeaf, fp);

  // ファイルクローズ
  fclose(fp);

  // エンディアン
  if( isNeedSwap )
  {
    for(std::vector<Pedigree>::iterator it = pedigrees.begin(); it != pedigrees.end(); ++it)
    {
      BSwap64(&(*it));
    }
  }

  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// 木情報ファイルのヘッダー読み込み
cpm_ErrorCode
cpm_VoxelInfoLMR::LoadOctreeHeader( std::string octFile, OctHeader &header )
{
  cpm_ErrorCode ret = CPM_SUCCESS;

  // ファイルオープン
  FILE *fp = fopen( octFile.c_str(), "rb" );
  if( !fp )
  {
    return CPM_ERROR_LMR_OPEN_OCTFILE;
  }

  // 読み込み
  bool isNeedSwap;
  ret = LoadOctreeHeader( fp, header, isNeedSwap );

  // ファイルクローズ
  fclose(fp);
  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// 木情報ファイルのヘッダー読み込み
cpm_ErrorCode
cpm_VoxelInfoLMR::LoadOctreeHeader( FILE *fp, OctHeader &header, bool &isNeedSwap )
{
  cpm_ErrorCode ret = CPM_SUCCESS;

  // ファイルポインタのチェック
  if( !fp )
  {
    return CPM_ERROR_LMR_OPEN_OCTFILE;
  }

  // ヘッダー読み込み
  isNeedSwap = false;
  fread(&header, sizeof(header), 1, fp);

  // エンディアン識別子のチェック
  if( header.identifier != OCTREE_FILE_IDENTIFIER )
  {
    BSwap32(&header.identifier);

    if( header.identifier != OCTREE_FILE_IDENTIFIER )
    {
      // スワップしても一致しなかった
      return CPM_ERROR_LMR_INVALID_OCTFILE;
    }

    // ここに来たときはエンディアン変換が必要
    isNeedSwap = true;
    for(int i = 0; i < 3; i++)
    {
      BSwap64(&header.org[i]);
      BSwap64(&header.rgn[i]);
      BSwap32(&header.rootDims[i]);
    }
    BSwap32(&header.maxLevel);
    BSwap64(&header.numLeaf);
  }

  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// グローバルの領域情報をセット
void
cpm_VoxelInfoLMR::SetGlobaliDomainInfo( S_OCT_DOMAIN_INFO &domainInfo )
{
  // 原点と領域サイズ
  m_globalDomainInfo.SetOrigin( domainInfo.origin );
  m_globalDomainInfo.SetRegion( domainInfo.region );

  // 1ルートあたりの自身のリーフレベルでの領域分割数
  int divNum = 1 << m_node->getLevel();

  // ピッチ、格子数、領域分割数(自身のリーフレベルにおける格子数、ピッチ、領域分割数とする)
  int vox[3], div[3];
  double pch[3];
  for( int i=0;i<3;i++ )
  {
    div[i] = divNum * m_octHeader.rootDims[i];
    vox[i] = div[i] * domainInfo.size[i];
    pch[i] = domainInfo.region[i] / double(vox[i]);
  }
  m_globalDomainInfo.SetVoxNum( vox );
  m_globalDomainInfo.SetPitch ( pch );
  m_globalDomainInfo.SetDivNum( div );

  // 活性サブドメイン情報はセットしない
  return;
}

////////////////////////////////////////////////////////////////////////////////
// ローカルの領域情報をセット
void
cpm_VoxelInfoLMR::SetLocalDomainInfo( S_OCT_DOMAIN_INFO &domainInfo )
{
  // 自身のぺディグリーを取得
  const Pedigree& pedigree = m_node->getPedigree();

  // 自身のレベルでの領域分割数(1ルートあたり)
  int div = pedigree.getUpperBound();

  // vox
  int vox[3];
  for( int m=0;m<3;m++ )
  {
    vox[m] = domainInfo.size[m];
  }
  m_localDomainInfo.SetVoxNum(vox);

  // origin
  Vec3d oct_origin = m_octree->getOrigin( m_node ); // 0.0～1.0の範囲でのリーフ原点座標
  double org[3];
  for( int m=0;m<3;m++ )
  {
    org[m] = oct_origin[m] * m_globalDomainInfo.GetRegion()[m] + m_globalDomainInfo.GetOrigin()[m];
  }
  m_localDomainInfo.SetOrigin(org);

  // region
  double rgn[3];
  for( int m=0;m<3;m++ )
  {
    rgn[m] = m_globalDomainInfo.GetRegion()[m] / double(m_octHeader.rootDims[m]) / double(div);
  }
  m_localDomainInfo.SetRegion(rgn);

  // pitch
  double pch[3];
  for( int m=0;m<3;m++ )
  {
    pch[m] = rgn[m] / double(vox[m]);
  }
  m_localDomainInfo.SetPitch(pch);

  // pos
  // 自身のリーフレベルでの分割位置をセット
  int pos[3];
  int divNum = 1 << m_node->getLevel();
  int rootID = pedigree.getRootID();
  const RootGrid *rootGrid = m_octree->getRootGrid();
  int ix = rootGrid->rootID2indexX(rootID);
  int iy = rootGrid->rootID2indexY(rootID);
  int iz = rootGrid->rootID2indexZ(rootID);
  int lx = pedigree.getX();
  int ly = pedigree.getY();
  int lz = pedigree.getZ();
  pos[0] = ix * divNum + lx;
  pos[1] = iy * divNum + ly;
  pos[2] = iz * divNum + lz;
  m_localDomainInfo.SetPos(pos);

  // head/tail
  // 自身のリーフレベルでのhead/Tailをセット
  for( int m=0;m<3;m++ )
  {
    m_voxelHeadIndex[m] = pos[m] * m_localDomainInfo.GetVoxNum()[m];
    m_voxelTailIndex[m] = m_voxelHeadIndex[m] + m_localDomainInfo.GetVoxNum()[m] - 1;
  }
}

////////////////////////////////////////////////////////////////////////////////
// 隣接情報の取得
void
cpm_VoxelInfoLMR::SetNeighborInfo()
{
  Partition part(m_nRank, m_octHeader.numLeaf);
  RootGrid *rootGrid = (RootGrid*)m_octree->getRootGrid();

  int cpm_face[6] = {X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS};
  int bcm_face[6] = {X_M    , X_P   , Y_M    , Y_P   , Z_M    , Z_P   };

  // 通常の隣接情報の生成
  {
    rootGrid->clearPeriodicX();
    rootGrid->clearPeriodicY();
    rootGrid->clearPeriodicZ();
    NeighborInfo *nInfo = m_octree->makeNeighborInfo( m_node, &part );
    for( int m=0;m<6;m++ )
    {
      // 隣接領域のレベル差
//      m_neighborLevelDiff[cpm_face[m]] = nInfo[bcm_face[m]].getLevelDifference();

      // 隣接ランク番号をセット
      int cnt=0;
      for( int i=0;i<4;i++ )
      {
        int nID = nInfo[bcm_face[m]].getID(Subface(i));
        if( nID < 0 ) continue;
        if( !nInfo[bcm_face[m]].isOuterBoundary() )
        {
          m_neighborRankID_LMR[cpm_face[m]][cnt++] = nID;
        }
      }
      m_neighborRankID[cpm_face[m]] = m_neighborRankID_LMR[cpm_face[m]][0];
    }
    delete [] nInfo;
  }

  // 周期境界の隣接情報の生成
  {
    rootGrid->setPeriodicX();
    rootGrid->setPeriodicY();
    rootGrid->setPeriodicZ();
    NeighborInfo *pInfo = m_octree->makeNeighborInfo( m_node, &part );
    for( int m=0;m<6;m++ )
    {
      // 隣接領域のレベル差
      m_neighborLevelDiff[cpm_face[m]] = pInfo[bcm_face[m]].getLevelDifference();

      // 周期境界で無い隣接が存在する場合、周期境界の隣接はセットしない
      if( !IsRankNull(m_neighborRankID[cpm_face[m]]) ) continue;

      // 隣接ランク番号をセット
      int cnt=0;
      for( int i=0;i<4;i++ )
      {
        int nID = pInfo[bcm_face[m]].getID(Subface(i));
        if( nID < 0 ) continue;
        m_periodicRankID_LMR[cpm_face[m]][cnt++] = nID;
      }
      m_periodicRankID[cpm_face[m]] = m_periodicRankID_LMR[cpm_face[m]][0];
    }
    delete [] pInfo;
  }

  // 周期境界フラグを元に戻しておく
  rootGrid->clearPeriodicX();
  rootGrid->clearPeriodicY();
  rootGrid->clearPeriodicY();
}

////////////////////////////////////////////////////////////////////////////////
// 木情報ファイルからリーフ数を取得する
int
cpm_VoxelInfoLMR::GetNumLeaf( std::string treeFile )
{
  // 領域情報の読み込み
  S_OCT_DOMAIN_INFO domainInfo;
  if( cpm_TextParserDomainLMR::Read(treeFile, domainInfo) != CPM_SUCCESS )
  {
    return 0;
  }

  // 木情報ファイルのヘッダーを読み込み
  OctHeader header;
  if( LoadOctreeHeader(domainInfo.octFile, header) != CPM_SUCCESS )
  {
    return 0;
  }

  // リーフ数
  return header.numLeaf;
}


////////////////////////////////////////////////////////////////////////////////
// 自ランクの隣接ランク番号を取得
const int*
cpm_VoxelInfoLMR::GetNeighborRankList( cpm_FaceFlag face, int &num ) const
{
  num = 0;
  for( int i=0;i<4;i++ )
  {
    if( !IsRankNull(m_neighborRankID_LMR[face][i]) ) num++;
  }
  return m_neighborRankID_LMR[face];
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの周期境界の隣接ランク番号を取得
const int*
cpm_VoxelInfoLMR::GetPeriodicRankList( cpm_FaceFlag face, int &num ) const
{
  num = 0;
  for( int i=0;i<4;i++ )
  {
    if( !IsRankNull(m_periodicRankID_LMR[face][i]) ) num++;
  }
  return m_periodicRankID_LMR[face];
}

////////////////////////////////////////////////////////////////////////////////
// 指定面におけるレベル差を取得
int
cpm_VoxelInfoLMR::GetNeighborLevelDiff( cpm_FaceFlag face ) const
{
  return m_neighborLevelDiff[face];
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの境界が外部境界かどうかを判定
bool
cpm_VoxelInfoLMR::IsOuterBoundary( cpm_FaceFlag face ) const
{
  // 隣のドメインが存在する
  int num = 0;
  const int *nlist = GetNeighborRankList( face, num );
  if( num > 0 )
  {
    return false;
  }

  // 分割数とポジション
  const int *div = GetDivNum();
  const int *pos = GetDivPos();

  // 分割数とポジションから外部境界かどうか判定
  if( face == X_MINUS && pos[0] == 0        ) return true;
  if( face == X_PLUS  && pos[0] == div[0]-1 ) return true;
  if( face == Y_MINUS && pos[1] == 0        ) return true;
  if( face == Y_PLUS  && pos[1] == div[1]-1 ) return true;
  if( face == Z_MINUS && pos[2] == 0        ) return true;
  if( face == Z_PLUS  && pos[2] == div[2]-1 ) return true;

  // ここに到達した場合は内部境界
  return false; 
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの境界が内部境界(隣が不活性ドメイン)かどうかを判定
bool
cpm_VoxelInfoLMR::IsInnerBoundary( cpm_FaceFlag face ) const
{
  // 隣のドメインが存在する
  int num = 0;
  const int *nlist = GetNeighborRankList( face, num );
  if( num > 0 )
  {
    return false;
  }

  // 分割数とポジション
  const int *div = GetDivNum();
  const int *pos = GetDivPos();

  // 分割数とポジションから外部境界かどうか判定
  if( face == X_MINUS && pos[0] == 0        ) return false;
  if( face == X_PLUS  && pos[0] == div[0]-1 ) return false;
  if( face == Y_MINUS && pos[1] == 0        ) return false;
  if( face == Y_PLUS  && pos[1] == div[1]-1 ) return false;
  if( face == Z_MINUS && pos[2] == 0        ) return false;
  if( face == Z_PLUS  && pos[2] == div[2]-1 ) return false;

  // ここに到達した場合は内部境界
  return true;
}

