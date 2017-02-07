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
 * @file   cpm_TextParserDomain.cpp
 * LMR用領域情報のTextParserクラスのソースファイル
 * @date   2015/03/27
 */
#include "cpm_TextParserDomainLMR.h"
#include "cpm_PathUtil.h"

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_TextParserDomainLMR::cpm_TextParserDomainLMR()
  : cpm_TextParser()
{
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_TextParserDomainLMR::~cpm_TextParserDomainLMR()
{
}

////////////////////////////////////////////////////////////////////////////////
// 読み込み(静的関数)
cpm_ErrorCode
cpm_TextParserDomainLMR::Read( std::string filename, S_OCT_DOMAIN_INFO &domainInfo )
{
  // インスタンス
  cpm_TextParserDomainLMR tp;

  // 読み込みメイン
  return tp.ReadMain( filename, domainInfo );
}

////////////////////////////////////////////////////////////////////////////////
// 読み込み処理のメイン
cpm_ErrorCode
cpm_TextParserDomainLMR::ReadMain( std::string filename, S_OCT_DOMAIN_INFO &domainInfo )
{
  int errorcode = TP_NO_ERROR;

  // デフォルトの読み込み処理
  if( (errorcode = cpm_TextParser::Read( filename )) != TP_NO_ERROR )
  {
    return CPM_ERROR_TP_LMR_DOMAINFILE;
  }

  // Domainの読み込み
  if( (errorcode = ReadDomain( domainInfo )) != TP_NO_ERROR )
  {
    return CPM_ERROR_TP_LMR_DOMAIN;
  }

  // BCMTreeの読み込み
  if( (errorcode = ReadBCMTree( domainInfo, filename )) != TP_NO_ERROR )
  {
    return CPM_ERROR_TP_LMR_BCMTREE;
  }

  // LeafBlockの読み込み
  if( (errorcode = ReadLeafBlock( domainInfo )) != TP_NO_ERROR )
  {
    return CPM_ERROR_TP_LMR_LEAFBLOCK;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Domainの読み込み
int
cpm_TextParserDomainLMR::ReadDomain( S_OCT_DOMAIN_INFO &domainInfo )
{
  int ret;

  // カレントノードを変更
  std::string oldpos;
  if( (ret = m_tp->currentNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }
  if( (ret = m_tp->changeNode( "/Domain" )) != TP_NO_ERROR )
  {
    return ret;
  }

  // リーフのラベルリストを取得
  std::vector<std::string> labels;
  if( (ret = m_tp->getLabels( labels )) != TP_NO_ERROR )
  {
    return ret;
  }

  bool borg = false;
  bool brgn = false;
  for( size_t i=0;i<labels.size();i++ )
  {
    std::string label = labels[i];

    // G_origin
    if( !borg && cpm_strCompare( label, "GlobalOrigin" ) == 0 )
    {
      if( readVector( label, domainInfo.origin, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_ORG;
      }
      borg = true;
      continue;
    }

    // G_region
    if( !brgn && cpm_strCompare( label, "GlobalRegion" ) == 0 )
    {
      if( readVector( label, domainInfo.region, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_RGN;
      }
      brgn = true;
      continue;
    }
  }

  // 取得チェック
  if( !borg )
  {
    return CPM_ERROR_TP_INVALID_G_ORG;
  }
  if( !brgn )
  {
    return CPM_ERROR_TP_INVALID_G_RGN;
  }

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
cpm_TextParserDomainLMR::ReadBCMTree( S_OCT_DOMAIN_INFO &domainInfo, std::string tpFile )
{
  int ret;

  // カレントノードを変更
  std::string oldpos;
  if( (ret = m_tp->currentNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }
  if( (ret = m_tp->changeNode( "BCMTree" )) != TP_NO_ERROR )
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
    if( !bfname && cpm_strCompare( label, "TreeFile" ) == 0 )
    {
      if( (ret = m_tp->getValue( label, filename )) != TP_NO_ERROR )
      {
        return ret;
      }
      bfname = true;
#ifdef _DEBUG
      std::cout << "TreeFile filename = " << filename << std::endl;
#endif
      break;
    }
  }

  // 取得チェック
  if( !bfname )
  {
    return CPM_ERROR_TP_LMR_BCMTREE;
  }

  // パスを含めたファイル名をセット
  if( bfname )
  {
    std::string fname = filename;
    // 相対パスのとき、元のtpファイルからの相対パスとする
    if( !CPM_PATH::cpmPath_isAbsolute( filename ) )
    {
      std::string dirName = CES::DirName(tpFile);
      fname = CPM_PATH::cpmPath_concat( dirName, fname );
#if 0
      std::cout << filename << " is Absolute path" << std::endl;
      std::cout << "dir name = " << dirName << std::endl;
      std::cout << "fname = " << fname << std::endl;
#endif
    }
    domainInfo.octFile = fname;
  }

  // 位置を元に戻す
  if( (ret = m_tp->changeNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// LeafBlockの読み込み
int
cpm_TextParserDomainLMR::ReadLeafBlock( S_OCT_DOMAIN_INFO &domainInfo )
{
  int ret;

  std::string oldpos;
  if( (ret = m_tp->currentNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }

  // LeafBlock直下
  {
    // カレントノードを変更
    if( (ret = m_tp->changeNode( "/LeafBlock" )) != TP_NO_ERROR )
    {
      return ret;
    }

    // リーフのラベルリストを取得
    std::vector<std::string> labels;
    if( (ret = m_tp->getLabels( labels )) != TP_NO_ERROR )
    {
      return ret;
    }

    bool bsize = false;
    for( size_t i=0;i<labels.size();i++ )
    {
      std::string label = labels[i];

      // Size
      if( !bsize && cpm_strCompare( label, "Size" ) == 0 )
      {
        if( readVector( label, domainInfo.size, 3 ) != TP_NO_ERROR )
        {
          return CPM_ERROR_TP_LMR_LEAFBLOCK;
        }
        bsize = true;
        break;
      }
    }

    // 取得チェック
    if( !bsize )
    {
      return CPM_ERROR_TP_LMR_LEAFBLOCK;
    }

    // 偶数チェック
    if( domainInfo.size[0]%2 != 0 ||
        domainInfo.size[1]%2 != 0 ||
        domainInfo.size[2]%2 != 0 )
    {
      return CPM_ERROR_TP_LMR_SIZE_NOT_EVEN;
    }
  }

  // Unit
  {
    // カレントノードを変更
    if( (ret = m_tp->changeNode( "/LeafBlock/Unit" )) != TP_NO_ERROR )
    {
      return ret;
    }

    // リーフのラベルリストを取得
    std::vector<std::string> labels;
    if( (ret = m_tp->getLabels( labels )) != TP_NO_ERROR )
    {
      return ret;
    }

    // 読み込み
    for( size_t i=0;i<labels.size();i++ )
    {
      std::string label = labels[i];

      // Length
      if( cpm_strCompare( label, "Length" ) == 0 )
      {
        if( (ret = m_tp->getValue( label, domainInfo.unitLength )) != TP_NO_ERROR )
        {
          return CPM_ERROR_TP_LMR_UNIT;
        }
        break;
      }
    }
  }

  // 位置を元に戻す
  if( (ret = m_tp->changeNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }

  return CPM_SUCCESS;
}
