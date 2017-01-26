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
 * @file   cpm_BaseParaManager_MPI.cpp
 * パラレルマネージャクラスのMPIインターフェイス関数ソースファイル
 * @date   2012/05/31
 */
#include "stdlib.h"
#include "cpm_BaseParaManager.h"

#if !defined(_WIN32) && !defined(WIN32)
#include <unistd.h> // for gethostname() of FX10/K
#endif

////////////////////////////////////////////////////////////////////////////////
// MPI_Datatypeを取得
MPI_Datatype
cpm_BaseParaManager::GetMPI_Datatype(int datatype)
{
  if( datatype == CPM_REAL )
  {
    if( RealIsDouble() ) return MPI_DOUBLE;
    return MPI_FLOAT;
  }
  else if( datatype == CPM_CHAR )               return MPI_CHAR;
  else if( datatype == CPM_SHORT )              return MPI_SHORT;
  else if( datatype == CPM_INT )                return MPI_INT;
  else if( datatype == CPM_LONG )               return MPI_LONG;
  else if( datatype == CPM_FLOAT )              return MPI_FLOAT;
  else if( datatype == CPM_DOUBLE )             return MPI_DOUBLE;
  else if( datatype == CPM_LONG_DOUBLE )        return MPI_LONG_DOUBLE;
  else if( datatype == CPM_UNSIGNED_CHAR )      return MPI_UNSIGNED_CHAR;
  else if( datatype == CPM_UNSIGNED_SHORT )     return MPI_UNSIGNED_SHORT;
  else if( datatype == CPM_UNSIGNED )           return MPI_UNSIGNED;
  else if( datatype == CPM_UNSIGNED_LONG )      return MPI_UNSIGNED_LONG;
#ifdef MPI_LONG_LONG_INT
  else if( datatype == CPM_LONG_LONG_INT )      return MPI_LONG_LONG_INT;
#endif
#ifdef MPI_LONG_LONG
  else if( datatype == CPM_LONG_LONG )          return MPI_LONG_LONG;
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( datatype == CPM_UNSIGNED_LONG_LONG ) return MPI_UNSIGNED_LONG_LONG;
#endif

  return MPI_DATATYPE_NULL;
}

////////////////////////////////////////////////////////////////////////////////
// MPI_Opを取得
MPI_Op
cpm_BaseParaManager::GetMPI_Op(int op)
{
  if     ( op == CPM_MAX    ) return MPI_MAX;
  else if( op == CPM_MIN    ) return MPI_MIN;
  else if( op == CPM_SUM    ) return MPI_SUM;
  else if( op == CPM_PROD   ) return MPI_PROD;
  else if( op == CPM_LAND   ) return MPI_LAND;
  else if( op == CPM_BAND   ) return MPI_BAND;
  else if( op == CPM_LOR    ) return MPI_LOR;
  else if( op == CPM_BOR    ) return MPI_BOR;
  else if( op == CPM_LXOR   ) return MPI_LXOR;
  else if( op == CPM_BXOR   ) return MPI_BXOR;
//  else if( op == CPM_MINLOC ) return CPM_MINLOC; // not support
//  else if( op == CPM_MAXLOC ) return CPM_MAXLOC; // not support

  return MPI_OP_NULL;
}

////////////////////////////////////////////////////////////////////////////////
// ランク番号の取得
int
cpm_BaseParaManager::GetMyRankID( int procGrpNo )
{
  // 不正なプロセスグループ番号
  if( procGrpNo < 0 || procGrpNo >= int(m_procGrpList.size()) )
  {
    // プロセスグループが存在しない
    return getRankNull();
  }

  // コミュニケータをチェック
  MPI_Comm comm = m_procGrpList[procGrpNo];
  if( IsCommNull(comm) )
  {
    // プロセスグループに自ランクが含まれない
   return getRankNull();
  }

  // ランク番号を取得
  int rankNo;
  MPI_Comm_rank( comm, &rankNo );

  // ランク番号
  return rankNo;
}

////////////////////////////////////////////////////////////////////////////////
// ランク数の取得
int
cpm_BaseParaManager::GetNumRank( int procGrpNo )
{
  // 不正なプロセスグループ番号
  if( procGrpNo < 0 || procGrpNo >= int(m_procGrpList.size()) )
  {
    // プロセスグループが存在しない
    return -1;
  }

  // コミュニケータをチェック
  MPI_Comm comm = m_procGrpList[procGrpNo];
  if( IsCommNull(comm) )
  {
    // プロセスグループに自ランクが含まれない
    return -1;
  }

  // ランク数を取得
  int nrank;
  MPI_Comm_size( comm, &nrank );

  // ランク数
  return nrank;
}

////////////////////////////////////////////////////////////////////////////////
// ホスト名の取得
std::string
cpm_BaseParaManager::GetHostName()
{
  char name[512];
  memset(name, 0x00, sizeof(char)*512);
  if( gethostname(name, 512) != 0 ) return std::string("");
  return std::string(name);
}

////////////////////////////////////////////////////////////////////////////////
// コミュニケータの取得
MPI_Comm
cpm_BaseParaManager::GetMPI_Comm( int procGrpNo )
{
  // 不正なプロセスグループ番号
  if( procGrpNo < 0 || procGrpNo >= int(m_procGrpList.size()) )
  {
    // プロセスグループが存在しない
    return getCommNull();
  }

  return m_procGrpList[procGrpNo];
}

////////////////////////////////////////////////////////////////////////////////
// Abort
void
cpm_BaseParaManager::Abort( int errorcode )
{
  // MPI_Abort
  MPI_Abort(MPI_COMM_WORLD, errorcode);
  exit(errorcode);
}

////////////////////////////////////////////////////////////////////////////////
// Barrier
cpm_ErrorCode
cpm_BaseParaManager::Barrier( int procGrpNo )
{
  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm( procGrpNo );
  if( IsCommNull( comm ) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Barrier
  if( MPI_Barrier(comm) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_BARRIER;
  }
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Wait
cpm_ErrorCode
cpm_BaseParaManager::Wait( MPI_Request *request )
{
  if( !request )
  {
    return CPM_ERROR_INVALID_PTR;
  }
  if( *request == MPI_REQUEST_NULL )
  {
    return CPM_ERROR_MPI_INVALID_REQUEST;
  }

  // MPI_Wait
  MPI_Status status;
  if( MPI_Wait( request, &status ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_WAIT;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Waitall
cpm_ErrorCode
cpm_BaseParaManager::Waitall( int count, MPI_Request requests[] )
{
  // status
  int cnt = 0;
  MPI_Status  *stat = new MPI_Status[count];
  MPI_Request *req  = new MPI_Request[count];
  for( int i=0;i<count;i++ )
  {
    if( requests[i] != MPI_REQUEST_NULL ) req[cnt++] = requests[i];
  }
  if( cnt == 0 )
  {
    delete [] stat;
    delete [] req;
    return CPM_SUCCESS;
  }

  // MPI_Waitall
  if( MPI_Waitall( cnt, req, stat ) != MPI_SUCCESS )
  {
    delete [] stat;
    delete [] req;
    return CPM_ERROR_MPI_WAITALL;
  }

  delete [] stat;
  delete [] req;
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Bcast(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Bcast( MPI_Datatype dtype, void *buf, int count, int root, int procGrpNo )
{
  if( !buf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Bcast
  if( MPI_Bcast( buf, count, dtype, root, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_BCAST;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Send(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Send( MPI_Datatype dtype, void *buf, int count, int dest, int procGrpNo )
{
  if( !buf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Send
  int tag = 1;
  if( MPI_Send( buf, count, dtype, dest, tag, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_SEND;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Recv(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Recv( MPI_Datatype dtype, void *buf, int count, int source, int procGrpNo ) 
{
  if( !buf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Send
  int tag = 1;
  MPI_Status status;
  if( MPI_Recv( buf, count, dtype, source, tag, comm, &status ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_SEND;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Isend(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Isend( MPI_Datatype dtype, void *buf, int count, int dest, MPI_Request *request, int procGrpNo )
{
  if( !buf || !request )
  {
    return CPM_ERROR_INVALID_PTR;
  }
  *request = MPI_REQUEST_NULL;

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Isend
  int tag = 1;
  if( MPI_Isend( buf, count, dtype, dest, tag, comm, request ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_ISEND;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Irecv(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Irecv( MPI_Datatype dtype, void *buf, int count, int source, MPI_Request *request, int procGrpNo )
{
  if( !buf || !request )
  {
    return CPM_ERROR_INVALID_PTR;
  }
  *request = MPI_REQUEST_NULL;

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Irecv
  int tag = 1;
  if( MPI_Irecv( buf, count, dtype, source, tag, comm, request ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_IRECV;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Allreduce(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Allreduce( MPI_Datatype dtype, void *sendbuf, void *recvbuf, int count, MPI_Op op, int procGrpNo )
{
  if( !sendbuf || !recvbuf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Allreduce
  if( MPI_Allreduce( sendbuf, recvbuf, count, dtype, op, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_ALLREDUCE;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Gather(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Gather( MPI_Datatype stype, void *sendbuf, int sendcnt
                       , MPI_Datatype rtype, void *recvbuf, int recvcnt
                       , int root, int procGrpNo )
{
  if( !sendbuf || !recvbuf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  //コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Gather
  if( MPI_Gather( sendbuf, sendcnt, stype, recvbuf, recvcnt, rtype, root, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_GATHER;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Allgather(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Allgather( MPI_Datatype stype, void *sendbuf, int sendcnt
                          , MPI_Datatype rtype, void *recvbuf, int recvcnt
                          , int procGrpNo )
{
  if( !sendbuf || !recvbuf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  //コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Allaather
  if( MPI_Allgather( sendbuf, sendcnt, stype, recvbuf, recvcnt, rtype, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_ALLGATHER;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Gatherv(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Gatherv( MPI_Datatype stype, void *sendbuf, int sendcnt
                        , MPI_Datatype rtype, void *recvbuf, int *recvcnts
                        , int *displs, int root, int procGrpNo )
{
  if( !sendbuf || !recvbuf || !recvcnts || !displs )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  //コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Gather
  if( MPI_Gatherv( sendbuf, sendcnt, stype, recvbuf, recvcnts, displs
                 , rtype, root, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_GATHERV;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Allgatherv(MPI_Datatype指定)
cpm_ErrorCode
cpm_BaseParaManager::Allgatherv( MPI_Datatype stype, void *sendbuf, int sendcnt
                           , MPI_Datatype rtype, void *recvbuf, int *recvcnts
                           , int *displs, int procGrpNo )
{
  if( !sendbuf || !recvbuf || !recvcnts || !displs )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Allaather
  if( MPI_Allgatherv( sendbuf, sendcnt, stype, recvbuf, recvcnts, displs, rtype, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_ALLGATHERV;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Wait
cpm_ErrorCode
cpm_BaseParaManager::cpm_Wait( int reqNo )
{
  // MPI_Request
  MPI_Request *req = m_reqList.Get(reqNo);
  if( !req )
  {
    return CPM_ERROR_INVALID_OBJKEY;
  }

  // MPI_Wait
  MPI_Status status;
  if( MPI_Wait( req, &status ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_WAIT;
  }

  // 削除
  return m_reqList.Delete(reqNo);
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Waitall
cpm_ErrorCode
cpm_BaseParaManager::cpm_Waitall( int count, int reqNoList[] )
{
  // MPI_Request
  MPI_Status* stat = new MPI_Status [count];
  MPI_Request* req = new MPI_Request[count];
  int cnt = 0;
  for( int i=0;i<count;i++ )
  {
    MPI_Request *r = m_reqList.Get( reqNoList[i] );
    if( r )
    {
      r[cnt++] = *req;
    }
  }
  if( cnt == 0 )
  {
    delete [] stat;
    delete [] req;
    return CPM_ERROR_INVALID_OBJKEY;
  }

  // MPI_Waitall
  if( MPI_Waitall( cnt, req, stat ) != MPI_SUCCESS )
  {
    delete [] stat;
    delete [] req;
    return CPM_ERROR_MPI_WAITALL;
  }

  // 削除
  for( int i=0;i<count;i++ )
  {
    m_reqList.Delete( reqNoList[i] );
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Isend
cpm_ErrorCode
cpm_BaseParaManager::cpm_Isend( void *buf, int count, int datatype, int dest, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request *req = m_reqList.Create();
  if( !req )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // Isend
  cpm_ErrorCode ret = Isend( dtype, buf, count, dest, req, procGrpNo );
  if( ret != MPI_SUCCESS )
  {
    delete req;
    return ret;
  }

  // MPI_Requestを登録
  if( (*reqNo = m_reqList.Add(req) ) < 0 )
  {
    delete req;
    return CPM_ERROR_REGIST_OBJKEY;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Irecv
cpm_ErrorCode
cpm_BaseParaManager::cpm_Irecv( void *buf, int count, int datatype, int source, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Irecv
  MPI_Request req;
  cpm_ErrorCode ret = Irecv( dtype, buf, count, source, &req, procGrpNo );
  if( ret != MPI_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを登録
  MPI_Request *r = m_reqList.Create();
  *r = req;
  if( (*reqNo = m_reqList.Add(r) ) < 0 )
  {
    delete r;
    return CPM_ERROR_REGIST_OBJKEY;
  }

  return CPM_SUCCESS;
}

