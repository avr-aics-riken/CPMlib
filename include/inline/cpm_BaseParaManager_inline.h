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
 * @file   cpm_BaseParaManager_inline.h
 * パラレルマネージャクラスのinline関数ヘッダーファイル
 * @date   2012/05/31
 */

#ifndef _CPM_BASEPARAMANAGER_INLINE_H_
#define _CPM_BASEPARAMANAGER_INLINE_H_

////////////////////////////////////////////////////////////////////////////////
// 配列の初期化処理
template<class T> CPM_INLINE
void
cpm_BaseParaManager::InitArray( T *array, size_t size )
{
  size_t sz = size * sizeof(T);
  memset( array, 0x00, sz );
}

////////////////////////////////////////////////////////////////////////////////
// 配列のコピー
template<class T> CPM_INLINE
void
cpm_BaseParaManager::CopyArray( T *source, T *dist, size_t size )
{
  size_t sz = size * sizeof(T);
  memcpy( dist, source, sz );
}

////////////////////////////////////////////////////////////////////////////////
// MPI_Datatypeを取得
template<class T> CPM_INLINE
MPI_Datatype
cpm_BaseParaManager::GetMPI_Datatype( T *ptr )
{
  if( !ptr )
  {
    return MPI_DATATYPE_NULL;
  }

  if     ( typeid(ptr) == typeid(char*)  )              return MPI_CHAR;
  else if( typeid(ptr) == typeid(short*) )              return MPI_SHORT;
  else if( typeid(ptr) == typeid(int*) )                return MPI_INT;
  else if( typeid(ptr) == typeid(long*) )               return MPI_LONG;
  else if( typeid(ptr) == typeid(float*) )              return MPI_FLOAT;
  else if( typeid(ptr) == typeid(double*) )             return MPI_DOUBLE;
  else if( typeid(ptr) == typeid(long double*) )        return MPI_LONG_DOUBLE;
  else if( typeid(ptr) == typeid(unsigned char*) )      return MPI_UNSIGNED_CHAR;
  else if( typeid(ptr) == typeid(unsigned short*) )     return MPI_UNSIGNED_SHORT;
  else if( typeid(ptr) == typeid(unsigned*) )           return MPI_UNSIGNED;
  else if( typeid(ptr) == typeid(unsigned int*) )       return MPI_UNSIGNED;
  else if( typeid(ptr) == typeid(unsigned long*) )      return MPI_UNSIGNED_LONG;
#ifdef MPI_LONG_LONG_INT
  else if( typeid(ptr) == typeid(long long int*) )      return MPI_LONG_LONG_INT;
#endif
#ifdef MPI_LONG_LONG
  else if( typeid(ptr) == typeid(long long*) )          return MPI_LONG_LONG;
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( typeid(ptr) == typeid(unsigned long long*) ) return MPI_UNSIGNED_LONG_LONG;
#endif

  return MPI_DATATYPE_NULL;
}

////////////////////////////////////////////////////////////////////////////////
// Broadcast
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Bcast( T *buf, int count, int root, int procGrpNo )
{
  // 型を取得
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype(buf);
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Bcast
  return Bcast( dtype, (void*)buf, count, root, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// Send
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Send( T *buf, int count, int dest, int procGrpNo )
{
  // 型を取得
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype(buf);
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Send
  return Send( dtype, (void*)buf, count, dest, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// Recv
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Recv( T *buf, int count, int source, int procGrpNo )
{
  // 型を取得
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype(buf);
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Recv
  return Recv( dtype, (void*)buf, count, source, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// Isend
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Isend( T *buf, int count, int dest, MPI_Request *request, int procGrpNo )
{
  // 型を取得
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype(buf);
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Isend
  return Isend( dtype, (void*)buf, count, dest, request, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// Irecv
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Irecv( T *buf, int count, int source, MPI_Request *request, int procGrpNo )
{
  // 型を取得
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype(buf);
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Irecv
  return Irecv( dtype, (void*)buf, count, source, request, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// Allreduce
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Allreduce( T *sendbuf, T *recvbuf, int count, MPI_Op op, int procGrpNo )
{
  // 型を取得
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype(sendbuf);
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Allreduce
  return Allreduce( dtype, sendbuf, recvbuf, count, op, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// Gather
template<class Ts, class Tr> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Gather( Ts *sendbuf, int sendcnt, Tr *recvbuf, int recvcnt, int root, int procGrpNo )
{
  // 型を取得
  MPI_Datatype stype = cpm_BaseParaManager::GetMPI_Datatype(sendbuf);
  if( stype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }
  MPI_Datatype rtype = cpm_BaseParaManager::GetMPI_Datatype(recvbuf);
  if( rtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Gather
  return Gather( stype, sendbuf, sendcnt, rtype, recvbuf, recvcnt, root, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// Allgather
template<class Ts, class Tr> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Allgather( Ts *sendbuf, int sendcnt, Tr *recvbuf, int recvcnt, int procGrpNo )
{
  // 型を取得
  MPI_Datatype stype = cpm_BaseParaManager::GetMPI_Datatype(sendbuf);
  if( stype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }
  MPI_Datatype rtype = cpm_BaseParaManager::GetMPI_Datatype(recvbuf);
  if( rtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Allgather
  return Allgather( stype, (void*)sendbuf, sendcnt, rtype, (void*)recvbuf, recvcnt, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// Gatherv
template<class Ts, class Tr> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Gatherv( Ts *sendbuf, int sendcnt, Tr *recvbuf, int *recvcnts, int *displs, int root, int procGrpNo )
{
  // 型を取得
  MPI_Datatype stype = cpm_BaseParaManager::GetMPI_Datatype(sendbuf);
  if( stype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }
  MPI_Datatype rtype = cpm_BaseParaManager::GetMPI_Datatype(recvbuf);
  if( rtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Gatherv
  return Gatherv( stype, (void*)sendbuf, sendcnt, rtype, (void*)recvbuf, recvcnts, displs, root, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// Allgatherv
template<class Ts, class Tr> CPM_INLINE
cpm_ErrorCode
cpm_BaseParaManager::Allgatherv( Ts *sendbuf, int sendcnt, Tr *recvbuf, int *recvcnts, int *displs, int procGrpNo )
{
  // 型を取得
  MPI_Datatype stype = cpm_BaseParaManager::GetMPI_Datatype(sendbuf);
  if( stype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }
  MPI_Datatype rtype = cpm_BaseParaManager::GetMPI_Datatype(recvbuf);
  if( rtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Allaather
  return Allgatherv( stype, (void*)sendbuf, sendcnt, rtype, (void*)recvbuf, recvcnts, displs, procGrpNo );
}

#endif /* _CPM_BASEPARAMANAGER_INLINE_H_ */
