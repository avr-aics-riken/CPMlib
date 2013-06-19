/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_TextParserDomain.cpp
 * CPM領域情報のTextParserクラスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "cpm_TextParserDomain.h"
#include "cpm_PathUtil.h"

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_TextParserDomain::cpm_TextParserDomain()
  : cpm_TextParser()
{
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_TextParserDomain::~cpm_TextParserDomain()
{
}

////////////////////////////////////////////////////////////////////////////////
// 読み込み(静的関数)
cpm_GlobalDomainInfo*
cpm_TextParserDomain::Read( std::string filename, int &errorcode )
{
  // インスタンス
  cpm_TextParserDomain tp;

  // 読み込みメイン
  return tp.ReadMain( filename, errorcode );
}

////////////////////////////////////////////////////////////////////////////////
// 読み込み処理のメイン
cpm_GlobalDomainInfo*
cpm_TextParserDomain::ReadMain( std::string filename, int &errorcode )
{
  errorcode = TP_NO_ERROR;

  // デフォルトの読み込み処理
  if( (errorcode = cpm_TextParser::Read( filename )) != TP_NO_ERROR )
  {
    return NULL;
  }

  // 領域情報のインスタンス
  cpm_GlobalDomainInfo *dInfo = new cpm_GlobalDomainInfo();
  if( !dInfo )
  {
    errorcode = CPM_ERROR_INVALID_PTR;
    return NULL;
  }

  // DomainInfoの読み込み
  if( (errorcode = ReadDomainInfo( dInfo )) != TP_NO_ERROR )
  {
    delete dInfo;
    return NULL;
  }

  // SubdomainInfoの読み込み
  if( (errorcode = ReadSubdomainInfo( dInfo, filename )) != TP_NO_ERROR )
  {
    delete dInfo;
    return NULL;
  }

  return dInfo;
}

////////////////////////////////////////////////////////////////////////////////
// DomainInfoの読み込み
int
cpm_TextParserDomain::ReadDomainInfo( cpm_GlobalDomainInfo* dInfo )
{
  int ret;

  if( !dInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // カレントノードを変更
  std::string oldpos;
  if( (ret = m_tp->currentNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }
  if( (ret = m_tp->changeNode( "/DomainInfo" )) != TP_NO_ERROR )
  {
    return ret;
  }

  // リーフのラベルリストを取得
  std::vector<std::string> labels;
  if( (ret = m_tp->getLabels( labels )) != TP_NO_ERROR )
  {
    return ret;
  }

  double org[3] = {0.0,0.0,0.0}; bool borg = false;
  int    vox[3] = {0,0,0};             bool bvox = false;
  double pch[3] = {0.0,0.0,0.0}; bool bpch = false;
  double rgn[3] = {0.0,0.0,0.0}; bool brgn = false;
  int       div[3] = {0,0,0};             bool bdiv = false;

  for( size_t i=0;i<labels.size();i++ )
  {
    std::string label = labels[i];

    // G_origin
    if( !borg && cpm_strCompare( label, "G_origin" ) == 0 )
    {
      if( readVector( label, org, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_ORG;
      }
      borg = true;
      continue;
    }

    // G_voxel
    if( !bvox && cpm_strCompare( label, "G_voxel" ) == 0 )
    {
      if( readVector( label, vox, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_VOXEL;
      }
      bvox=true;
      continue;
    }

    // G_pitch
    if( !bpch && cpm_strCompare( label, "G_pitch" ) == 0 )
    {
      if( readVector( label, pch, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_PITCH;
      }
      bpch = true;
      continue;
    }

    // G_region
    if( !brgn && cpm_strCompare( label, "G_region" ) == 0 )
    {
      if( readVector( label, rgn, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_PITCH;
      }
      brgn = true;
      continue;
    }

    // G_div
    if( !bdiv && cpm_strCompare( label, "G_div" ) == 0 )
    {
      if( readVector( label, div, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_DIV;
      }
      bdiv = true;
      continue;
    }
  }

  // G_orgをセット(必須)
  if( borg )
    dInfo->SetOrigin( org );
  else
    return CPM_ERROR_TP_INVALID_G_ORG;

  // G_voxelをセット(オプション)
  if( bvox )
  {
    if( vox[0] <= 0 || vox[1] <= 0 || vox[2] <= 0 )
      return CPM_ERROR_TP_INVALID_G_VOXEL;
  }
  dInfo->SetVoxNum( vox );

  // G_region、G_pitchをセット(G_region優先)
  // regionが必ず決定するように指定されている必要がある
  // 1.G_regionが指定されているとき
  //   G_voxが指定されていればpitchを自動計算
  // 2.G_regionが指定されておらずG_pitchが指定されているとき
  //   G_voxが指定されていればregionを自動計算
  //   G_voxが指定されていないときはエラー
  if( brgn )
  {
    if( rgn[0] <= 0.0 || rgn[0] <= 0.0 || rgn[0] <= 0.0 )
    {
      return CPM_ERROR_TP_INVALID_G_RGN;
    }
    dInfo->SetRegion( rgn );
    if( bvox )
    {
      for( int i=0;i<3;i++ ) pch[i] = rgn[i] / double(vox[i]);
      dInfo->SetPitch( pch );
    }
  }
  else if( bpch && bvox)
  {
    if( pch[0] <= 0.0 || pch[0] <= 0.0 || pch[0] <= 0.0 )
    {
      return CPM_ERROR_TP_INVALID_G_PITCH;
    }
    for( int i=0;i<3;i++ ) rgn[i] = pch[i] * double(vox[i]);
    dInfo->SetPitch( pch );
    dInfo->SetRegion( rgn );
  }
  else
  {
    return CPM_ERROR_TP_INVALID_G_RGN;
  }

  // G_divをセット(オプション)
  if( bdiv )
  {
    if( div[0] <= 0 || div[1] <= 0 || div[2] <= 0 )
    return CPM_ERROR_TP_INVALID_G_DIV;
  }
  dInfo->SetDivNum( div );

  // 位置を元に戻す
  if( (ret = m_tp->changeNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// SubdomainInfoの読み込み
int
cpm_TextParserDomain::ReadSubdomainInfo( cpm_GlobalDomainInfo* dInfo, std::string tpfname )
{
  int ret;

  if( !dInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 領域分割数の取得
  const int *div = dInfo->GetDivNum();
  if( !div )
  {
    return CPM_ERROR_TP_INVALID_G_DIV;
  }

  // カレントノードを変更
  std::string oldpos;
  if( (ret = m_tp->currentNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }
  if( (ret = m_tp->changeNode( "SubdomainInfo" )) != TP_NO_ERROR )
  {
    return ret;
  }

  // リーフのラベルリストを取得
  std::vector<std::string> labels;
  if( (ret = m_tp->getLabels( labels )) != TP_NO_ERROR )
  {
    return ret;
  }

  // 各ラベルをパース
  std::string filename = ""; bool bfname = false;
  for( size_t i=0;i<labels.size();i++ )
  {
    // ラベル
    std::string label = labels[i];

    // filename
    if( !bfname && cpm_strCompare( label, "filename" ) == 0 )
    {
      if( (ret = m_tp->getValue( label, filename )) != TP_NO_ERROR )
      {
        return ret;
      }
      bfname = true;
#ifdef _DEBUG
      std::cout << "ActiveSubdomainInfo filename = " << filename << std::endl;
#endif
      break;
    }
  }

  // ActiveSubdomainファイルの読み込み
  if( bfname )
  {
    std::string fname = filename;
    // 相対パスのとき、元のtpファイルからの相対パスとする
    if( !CPM_PATH::cpmPath_isAbsolute( filename ) )
    {
      std::string dirName = CES::DirName(tpfname);
      fname = CPM_PATH::cpmPath_concat( dirName, fname );
#if 0
      std::cout << filename << " is Absolute path" << std::endl;
      std::cout << "dir name = " << dirName << std::endl;
      std::cout << "fname = " << fname << std::endl;
#endif
    }
    if( (ret = dInfo->ReadActiveSubdomainFile(fname)) != CPM_SUCCESS )
    {
      return ret;
    }
  }

  // 位置を元に戻す
  if( (ret = m_tp->changeNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }

  return CPM_SUCCESS;
}

