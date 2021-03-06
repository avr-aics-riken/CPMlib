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
 * @file   cpm_DomainInfo.h
 * 領域情報クラスのヘッダーファイル
 * @date   2012/05/31
 */

#ifndef _CPM_DomainInfo_H_
#define _CPM_DomainInfo_H_

#include <vector>
#include "cpm_Base.h"
#include "cpm_EndianUtil.h"

/** CPMの領域情報クラス */
class cpm_DomainInfo : public cpm_Base
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** コンストラクタ */
  cpm_DomainInfo();

  /** デストラクタ */
  virtual ~cpm_DomainInfo();

  /** 情報のクリア */
  virtual void clear();

  /** 原点のセット
   *  @param[in] org 原点情報
   */
  void SetOrigin( double org[3] );

  /** 原点の取得
   *  @return 原点情報実数配列のポインタ
   */
  const double* GetOrigin() const;

  /** ピッチのセット
   *  @param[in] pch ピッチ情報
   */
  void SetPitch( double pch[3] );

  /** ピッチの取得
   *  @return ピッチ情報実数配列のポインタ
   */
  const double* GetPitch() const;

  /** 空間サイズのセット
   *  @param[in] rgn 空間サイズ情報
   */
  void SetRegion( double rgn[3] );

  /** 空間サイズの取得
   *  @return 空間サイズ情報実数配列のポインタ
   */
  const double* GetRegion() const;

  /** VOXEL数のセット
   *  @param[in] vox VOXEL数情報
   */
  void SetVoxNum( int vox[3] );

  /** VOXEL数の取得
   *  @return VOXEL数情報実数配列のポインタ
   */
  const int* GetVoxNum() const;

// 2016/01/22 FEAST add.s
  /** 頂点数のセット
   *  @param[in] nod 頂点数情報
   */
  void SetNodNum( int nod[3] );

  /** 頂点数の取得
   *  @return 頂点数情報整数配列のポインタ
   */
  const int* GetNodNum() const;
// 2016/01/22 FEAST add.e

  /** 領域情報のチェック \n
   *  - VoxelInitを実行する上で必要な情報がセットされているかをチェックする。
   *
   *  @return   終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode CheckData();

private:


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:
  double m_origin[3]; ///< 原点
  double m_region[3]; ///< 空間サイズ
  double m_pitch[3];  ///< ピッチ
  int    m_voxNum[3]; ///< VOXEL数

// 2016/01/22 FEAST add.s
  int    m_nodNum[3]; ///< 頂点数
// 2016/01/22 FESAT add.e

};

/** CPMのサブ領域情報クラス */
class cpm_ActiveSubdomainInfo : public cpm_Base
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** デフォルトコンストラクタ */
  cpm_ActiveSubdomainInfo();

  /** コンストラクタ
   *  @param[in] pos  領域分割内での位置
   */
  cpm_ActiveSubdomainInfo( int pos[3] );

  /** デストラクタ */
  virtual ~cpm_ActiveSubdomainInfo();

  /** 情報のクリア */
  virtual void clear();

  /** 位置のセット
   *  @param[in] pos  領域分割内での位置
   */
  void SetPos( int pos[3] );

  /** 位置の取得
   *  @return 位置情報整数配列のポインタ
   */
  const int* GetPos() const;

  /** 比較演算子
   *  @param[in] dom   比較対象の活性サブドメイン情報
   *  @retval    true  同じ位置情報を持つ
   *  @retval    false 違う位置情報を持つ
   */
  bool operator==( cpm_ActiveSubdomainInfo dom );

  /** 比較演算子
   *  @param[in] dom   比較対象の活性サブドメイン情報
   *  @retval    true  違う位置情報を持つ
   *  @retval    false 同じ位置情報を持つ
   */
  bool operator!=( cpm_ActiveSubdomainInfo dom );


private:


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:
  int m_pos[3];  ///< 領域分割内での位置


};

/** CPMの全体領域情報クラス */
class cpm_GlobalDomainInfo : public cpm_DomainInfo
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** コンストラクタ */
  cpm_GlobalDomainInfo();

  /** デストラクタ */
  virtual ~cpm_GlobalDomainInfo();

  /** 情報のクリア */
  virtual void clear();

  /** 領域分割数のセット
   *  @param[in] div 領域分割数
   */
  void SetDivNum( int div[3] );

  /** 領域分割数の取得
   *  @return 領域分割数整数配列のポインタ
   */
  const int* GetDivNum() const;

  /** 活性サブドメイン情報の存在チェック
   *  @param[in] subDomain チェックする活性サブドメイン情報
   *  @retval    true      存在する
   *  @retval    false     存在しない
   */
  bool IsExistSubdomain( cpm_ActiveSubdomainInfo subDomain );

  /** 活性サブドメイン情報の追加
   *  @param[in] subDomain 追加する活性サブドメイン情報
   *  @retval    true      追加した
   *  @retval    false     追加に失敗(同じ領域分割位置で追加済み)
   */
  bool AddSubdomain( cpm_ActiveSubdomainInfo subDomain );

  /** 活性サブドメインの数を取得
   *  - 活性サブドメインの数＝活性サブドメイン情報配列のサイズだが、
   *  この配列が空のとき、領域分割数でサブドメイン数を決定して返す
   *
   *  @return 活性サブドメインの数
   */
  int GetSubdomainNum() const;

  /** 活性サブドメインの数を取得(情報数)
   *  - 活性サブドメインの数＝活性サブドメイン情報配列のサイズ
   *
   *  @return 活性サブドメイン情報配列サイズ
   */
  int GetSubdomainArraySize() const;

  /** 活性サブドメイン情報を取得
   *  @param[in]  idx 登録順番号
   *  @return 活性サブドメイン情報ポインタ
   */
  const cpm_ActiveSubdomainInfo* GetSubdomainInfo( size_t idx ) const;

  /** 領域情報のチェック
   *  - VoxelInitを実行する上で必要な情報がセットされているかをチェックする。
   *   活性サブドメイン配列が空のとき、全領域が活性サブドメインになるため、
   *  このチェック関数内で活性サブドメイン情報を生成する。
   *
   *  @param[in] nRank 並列プロセス数
   *  @return   終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode CheckData( int nRank );

  /** ActiveSubdomainファイルの読み込み
   *  - ActiveSubdomainファイルを読み込み、活性ドメイン情報を生成する
   *
   *  @param[in]  subDomainFile ActiveSubdomainファイル名
   *  @return   終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode ReadActiveSubdomainFile( std::string subDomainFile );

  /** ActiveSubdomainファイルの読み込み(static関数)
   *  - ActiveSubdomainファイルを読み込み、活性ドメイン情報を生成する
   *
   *  @param[in]  subDomainFile ActiveSubdomainファイル名
   *  @param[out] subDomainInfo 活性ドメイン情報
   *  @param[out] div           ActiveSubdiomainファイル中の領域分割数
   *  @return   終了コード(CPM_SUCCESS=正常終了)
   */
  static cpm_ErrorCode ReadActiveSubdomainFile( std::string subDomainFile
                                              , std::vector<cpm_ActiveSubdomainInfo>& subDomainInfo
                                              , int div[3] );

  /** ActiveSubdomainファイルのエンディアンチェック
   *  - ActiveSubdomainファイルのエンディアンをチェック
   *
   *  @param[in]  ident               ActiveSubdomainファイルのIdentifier
   *  @retval     CPM_ENDIAN::Match   一致
   *  @retval     CPM_ENDIAN::UnMatch 不一致
   *  @retval     CPM_ENDIAN::UnKnown フォーマットが異なる
   */
  static CPM_ENDIAN::EMatchType isMatchEndianSbdmMagick( int ident );

private:


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:
  int m_divNum[3]; ///< 領域分割数
  std::vector<cpm_ActiveSubdomainInfo> m_subDomainInfo; ///< 活性サブドメイン情報


};

/** CPMのローカル領域情報クラス */
class cpm_LocalDomainInfo : public cpm_DomainInfo, public cpm_ActiveSubdomainInfo
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** コンストラクタ */
  cpm_LocalDomainInfo();

  /** デストラクタ */
  virtual ~cpm_LocalDomainInfo();

  /** 情報のクリア */
  virtual void clear();


private:


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:


};

#endif /* _CPM_DomainInfo_H_ */
