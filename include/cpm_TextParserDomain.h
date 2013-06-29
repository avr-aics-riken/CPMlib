/*
 * CPMlib - Cartesian Partition Manager Library
 *
 * Copyright (C) 2012-2013 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 */

/**
 * @file   cpm_TextParserDomain.h
 * 領域情報のテキストパーサークラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_TEXTPARSER_DOMAIN_H_
#define _CPM_TEXTPARSER_DOMAIN_H_

#include "cpm_TextParser.h"
#include "cpm_DomainInfo.h"

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
   *  @param[in]  filename 読み込むファイル名
   *  @param[out] errorcode CPMエラーコード
   *  @return 領域情報ポインタ
   */
  static cpm_GlobalDomainInfo* Read( std::string filename, int &errorcode ); 

  /** G_voxel,G_region,G_pitchからサイズを確定する
   *  本関数ではG_voxelを必須とし、G_regionとG_pitchの
   *  両方が指定されている場合、G_regionを優先する
   *  ユーザがこの基準を変更したい場合は、cpm_TextParserDomainクラスの
   *  派生クラスを作成し、本関数をオーバーライドし、その関数内に
   *  サイズ決定処理を記述する．
   *  その場合、Read関数もオーバーライドする必要がある．
   *  Example/cxx/user_TPdomain.h，CのTPdomainクラスを参照のこと．
   *  @param[in] bvox 領域分割情報ファイルのG_voxel記述フラグ(true:記述されている)
   *  @param[in] vox  領域分割情報ファイルのG_voxel(bvox=falseのとき不定値)
   *  @param[in] brgn 領域分割情報ファイルのG_region記述フラグ(true:記述されている)
   *  @param[in] rgn  領域分割情報ファイルのG_region(brgn=falseのとき不定値)
   *  @param[in] bpch 領域分割情報ファイルのG_pitch記述フラグ(true:記述されている)
   *  @param[in] pch  領域分割情報ファイルのG_pitch(bpch=falseのとき不定値)
   *  @return CPMエラーコード
   */
  virtual int DecideGregion( bool bvox, int    vox[3]
                            , bool brgn, double rgn[3]
                            , bool bpch, double pch[3] );
  
  
protected:
  /** 読み込み処理のメイン
   *  - TextParserクラスを用いて領域分割情報ファイルを読み込む
   *  - TextParserクラスのインスタンスはクリア(remove)される
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

