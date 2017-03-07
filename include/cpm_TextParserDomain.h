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
 * @file   cpm_TextParserDomain.h
 * 領域情報のテキストパーサークラスのヘッダーファイル
 * @date   2012/05/31
 */

#ifndef _CPM_TEXTPARSER_DOMAIN_H_
#define _CPM_TEXTPARSER_DOMAIN_H_

#include "cpm_TextParser.h"
#include "cpm_DomainInfo.h"
#include <string.h> // for strdup

/** CPMの領域情報テキストパーサークラス */
class cpm_TextParserDomain : public cpm_TextParser
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** コンストラクタ */
  cpm_TextParserDomain();

  /** デストラクタ */
  virtual ~cpm_TextParserDomain();

  /** 読み込み処理
   *  - TextParserクラスを用いて領域分割情報ファイルを読み込む
   *  - TextParserクラスのインスタンスはクリア(remove)される
   *
   *  @param[in]  filename 読み込むファイル名
   *  @param[out] errorcode CPMエラーコード
   *  @return 領域情報ポインタ
   */
  static cpm_GlobalDomainInfo* Read( std::string filename, int &errorcode );


private:
  /** 読み込み処理のメイン
   *  - TextParserクラスを用いて領域分割情報ファイルを読み込む
   *  - TextParserクラスのインスタンスはクリア(remove)される
   *
   *  @param[in]  filename 読み込むファイル名
   *  @param[out] errorcode CPMエラーコード
   *  @return 領域情報ポインタ
   */
  cpm_GlobalDomainInfo* ReadMain( std::string filename, int &errorcode );

  /** DomainInfoの読み込み
   *  @param[inout] dInfo 領域情報
   *  @return CPMエラーコード
   */
  int ReadDomainInfo( cpm_GlobalDomainInfo* dInfo );

  /** ActiveSubdomainInfoの読み込み
   *  @param[inout] dInfo   領域情報
   *  @param[in]    tpfname メインの領域分割情報ファイル名
   *  @return CPMエラーコード
   */
  int ReadSubdomainInfo( cpm_GlobalDomainInfo* dInfo, std::string tpfname );


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:


};

#endif /* _CPM_TEXTPARSER_DOMAIN_H_ */
