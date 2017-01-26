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
 * @file   cpm_TextParserDomainLMR.h
 * LMR用領域情報のテキストパーサークラスのヘッダーファイル
 * @date   2012/05/31
 */

#ifndef _CPM_TEXTPARSER_DOMAIN_LMR_H_
#define _CPM_TEXTPARSER_DOMAIN_LMR_H_

#include "cpm_TextParser.h"
#include <string.h> // for strdup

/** 領域情報ファイル構造体 */
struct S_OCT_DOMAIN_INFO
{
  double      origin[3];  ///< 原点座標
  double      region[3];  ///< 領域幅
  std::string octFile;    ///< octファイル名
  int         size[3];    ///< 1リーフの格子数
  std::string unitLength; ///< 長さ単位文字列

  /** コンストラクタ */
  S_OCT_DOMAIN_INFO()
  {
    origin[0] = origin[1] = origin[2] = 0.0;
    region[0] = region[1] = region[2] = 0.0;
    octFile = "";
    size[0] = size[1] = size[2] = 0.0;
    unitLength = "";
  }

  void print()
  {
    std::cout << "*** domain file info" << std::endl;
    std::cout << "origin  : " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
    std::cout << "region  : " << region[0] << " " << region[1] << " " << region[2] << std::endl;
    std::cout << "octfile : " << octFile << std::endl;
    std::cout << "size    : " << size[0] << " " << size[1] << " " << size[2] << std::endl;
    std::cout << "length  : " << unitLength << std::endl;
  }
};

/** LMR用領域情報テキストパーサークラス */
class cpm_TextParserDomainLMR : public cpm_TextParser
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** コンストラクタ */
  cpm_TextParserDomainLMR();

  /** デストラクタ */
  virtual ~cpm_TextParserDomainLMR();

  /** 読み込み処理
   *  - TextParserクラスを用いて領域分割情報ファイルを読み込む
   *  - TextParserクラスのインスタンスはクリア(remove)される
   *
   *  @param[in]  filename   読み込むファイル名
   *  @param[out] domainInfo 領域情報
   *  @return     CPMエラーコード
   */
  static cpm_ErrorCode Read( std::string filename, S_OCT_DOMAIN_INFO &domainInfo );


private:
  /** 読み込み処理のメイン
   *  - TextParserクラスを用いて領域分割情報ファイルを読み込む
   *  - TextParserクラスのインスタンスはクリア(remove)される
   *
   *  @param[in]    filename 読み込むファイル名
   *  @param[inout] domainInfo 領域情報
   *  @return       CPMエラーコード
   */
  cpm_ErrorCode ReadMain( std::string filename, S_OCT_DOMAIN_INFO &domainInfo ); 

  /** Domainの読み込み
   *  @param[inout] domainInfo 領域情報
   *  @return       CPMエラーコード
   */
  int ReadDomain( S_OCT_DOMAIN_INFO &domainInfo );

  /** BCMTreeの読み込み
   *  @param[inout] domainInfo 領域情報
   *  @param[in]    tpFile     元のテキストパーサー形式ファイル名
   *  @return       CPMエラーコード
   */
  int ReadBCMTree( S_OCT_DOMAIN_INFO &domainInfo, std::string tpFile );

  /** LeafBlockの読み込み
   *  @param[inout] domainInfo 領域情報
   *  @return       CPMエラーコード
   */
  int ReadLeafBlock( S_OCT_DOMAIN_INFO &domainInfo );


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:


};

#endif /* _CPM_TEXTPARSER_DOMAIN_LMR_H_ */

