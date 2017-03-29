#pragma once
#ifdef DISABLE_MPI
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <vector>
#include <time.h>

extern "C"
{

const int MPI_PROC_NULL		= -2;	///< null process
const int MPI_COMM_NULL		= -2;	///< null coommunicator
const int MPI_SUCCESS		= 0;	///< success code
const int MPI_COMM_WORLD	= 0;	///< world communicator
const int MPI_REQUEST_NULL	= 0;	///< null request

typedef int MPI_Comm;		///< mpi communicator
typedef int MPI_Group;		///< mpi group
typedef int MPI_Request;	///< mpi request
typedef int MPI_Status;		///< mpi status

static int cpm_MPIInitialized = 0;	///< initialized flag
static int cpm_MPIFinalized   = 0;	///< finalized flag

/// MPI data types
enum MPI_Datatype
{
  MPI_DATATYPE_NULL      =  0 ///< null
, MPI_CHAR               =  1 ///< char
, MPI_UNSIGNED_CHAR      =  2 ///< unsigned char
, MPI_BYTE               =  3 ///< byte(not support)
, MPI_SHORT              =  4 ///< short
, MPI_UNSIGNED_SHORT     =  5 ///< unsigned short
, MPI_INT                =  6 ///< int
, MPI_UNSIGNED           =  7 ///< unsigned
, MPI_LONG               =  8 ///< long
, MPI_UNSIGNED_LONG      =  9 ///< unsigned long
, MPI_FLOAT              = 10 ///< float
, MPI_DOUBLE             = 11 ///< double
, MPI_LONG_DOUBLE        = 12 ///< long double
, MPI_LONG_LONG_INT      = 13 ///< long long int
, MPI_LONG_LONG          = 13 ///< long long
, MPI_UNSIGNED_LONG_LONG = 51 ///< unsigned long long
, MPI_REAL               = 52 ///< REAL_TYPE
};

/// MPI opertor
enum MPI_Op
{
  MPI_OP_NULL = 0   ///< null
, MPI_MAX     = 100 ///< 最大値
, MPI_MIN     = 101 ///< 最小値
, MPI_SUM     = 102 ///< 和
, MPI_PROD    = 103 ///< 積
, MPI_LAND    = 104 ///< 論理積
, MPI_BAND    = 105 ///< ビット演算の積
, MPI_LOR     = 106 ///< 論理和
, MPI_BOR     = 107 ///< ビット演算の和
, MPI_LXOR    = 108 ///< 排他的論理和
, MPI_BXOR    = 109 ///< ビット演算の排他的論理和
, MPI_MINLOC  = 110 ///< 最大値と位置(not support)
, MPI_MAXLOC  = 111 ///< 最小値と位置(not support)
};

/// struct of send/recv information
struct CPM_STUBCOMMBUF_INFO
{
  int tag;		///< tag of mpi function argument
  MPI_Request req;	///< request for nonblocking
  size_t size;		///< size of data type
  size_t nmemb;		///< number of data
  void*  sbuf;		///< send data array
  void*  rbuf;		///< recv data array

  /// constractor
  CPM_STUBCOMMBUF_INFO()
  {
    tag = -1;
    req = MPI_REQUEST_NULL;
    size = 0;
    nmemb = 0;
    sbuf = NULL;
    rbuf = NULL;
  }

  /// destractor
  ~CPM_STUBCOMMBUF_INFO()
  {
    tag = -1;
    req = MPI_REQUEST_NULL;
    size = 0;
    nmemb = 0;
    sbuf = NULL;
    rbuf = NULL;
  }
};

/// send/recv buffer list (FIFO)
static std::vector<CPM_STUBCOMMBUF_INFO*> cpm_StubCommBufInfo;

/// get new request
static int cpm_StubGetRequest()
{
  int req = 0;
  bool exist = true;
  while(exist)
  {
    req++;
    exist = false;
    for( size_t i=0;i<cpm_StubCommBufInfo.size();i++ )
    {
      if( cpm_StubCommBufInfo[i]->req == req )
      {
        exist = true;
        break;
      }
    }
  }

  return req;
}

/// get size of MPI_Datatype(return value : byte)
static size_t cpm_StubGetDatatypeSize(MPI_Datatype type)
{
  if( type == MPI_CHAR )
    return sizeof(char);
  else if( type == MPI_UNSIGNED_CHAR )
    return sizeof(unsigned char);
  else if( type == MPI_BYTE )
    return sizeof(char);
  else if( type == MPI_SHORT )
    return sizeof(short);
  else if( type == MPI_UNSIGNED_SHORT )
    return sizeof(unsigned short);
  else if( type == MPI_INT )
    return sizeof(int);
  else if( type == MPI_UNSIGNED )
    return sizeof(unsigned);
  else if( type == MPI_LONG )
    return sizeof(long);
  else if( type == MPI_UNSIGNED_LONG )
    return sizeof(unsigned long);
  else if( type == MPI_FLOAT )
    return sizeof(float);
  else if( type == MPI_DOUBLE )
    return sizeof(double);
  else if( type == MPI_LONG_DOUBLE )
    return sizeof(long double);
  else if( type == MPI_LONG_LONG_INT )
    return sizeof(long long int);
  else if( type == MPI_LONG_LONG )
    return sizeof(long long);
  else if( type == MPI_UNSIGNED_LONG_LONG )
    return sizeof(unsigned long long);
  else if( type == MPI_REAL )
#ifdef _REAL_IS_DOUBLE_
    return sizeof(double);
#else
    return sizeof(float);
#endif

  return 1;
}

/// Returns an elapsed time on the calling processor 
static double MPI_Wtime()
{
  return clock();
}

/// Initialize the MPI execution environment 
static int MPI_Init(int *argc, char ***argv)
{
  cpm_MPIInitialized = 1;
  cpm_StubCommBufInfo.clear();
  return MPI_SUCCESS;
}

/// Indicates whether MPI_Init has been called. 
static int MPI_Initialized(int *flag)
{
  *flag = cpm_MPIInitialized;
  return MPI_SUCCESS;
}

/// Terminates MPI execution environment 
static int MPI_Finalize()
{
  cpm_MPIInitialized = 0;
  cpm_MPIFinalized = 1;
  for( size_t i=0;i<cpm_StubCommBufInfo.size();i++ )
  {
    delete cpm_StubCommBufInfo[i];
  }
  cpm_StubCommBufInfo.clear();
  return MPI_SUCCESS;
}

/// Indicates whether MPI_Finalize has been called. 
static int MPI_Finalized(int *flag)
{
  *flag = cpm_MPIFinalized;
  return MPI_SUCCESS;
}

/// Determines the size of the group associated with a communicator 
static int MPI_Comm_size(MPI_Comm comm, int *size)
{
  *size = 1;
  return MPI_SUCCESS;
}

/// Determines the rank of the calling process in the communicator 
static int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
  *rank = 0;
  return MPI_SUCCESS;
}

/// Accesses the group associated with given communicator 
static int MPI_Comm_group(MPI_Comm comm, MPI_Group *group)
{
  *group = 0;
  return MPI_SUCCESS;
}

/// Produces a group by reordering an existing group and taking only listed members 
static int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
{
  *newgroup = 0;
  return MPI_SUCCESS;
}

/// Creates a new communicator 
static int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)
{
  *newcomm = 0;
  return MPI_SUCCESS;
}

/// Terminates MPI execution environment 
static int MPI_Abort(MPI_Comm comm, int errorcode)
{
  exit(errorcode);
  return MPI_SUCCESS;
}

/// Blocks until all processes in the communicator have reached this routine. 
static int MPI_Barrier( MPI_Comm comm )
{
  return MPI_SUCCESS;
}

/// Broadcasts a message from the process with rank "root" to all other processes of the communicator 
static int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, 
               MPI_Comm comm )
{
  return MPI_SUCCESS;
}

/// Performs a blocking send 
static int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm)
{
  // 逐次なので自身(ランク0)のみ
  int myrank;
  MPI_Comm_rank(comm, &myrank);
  if( dest != myrank )
  {
    return MPI_SUCCESS;
  }

  // 同期通信のため、Sendが呼ばれたタイミングで受信情報が存在するはず

  // バッファを検索
  std::vector<CPM_STUBCOMMBUF_INFO*>::iterator it = cpm_StubCommBufInfo.begin();
  for( ;it!=cpm_StubCommBufInfo.end();it++ )
  {
    // tagが同じ受信情報を取得:
    CPM_STUBCOMMBUF_INFO *info = *it;
    if( info->tag == tag && info->rbuf )
    {
      // コピー
      size_t sz = cpm_StubGetDatatypeSize(datatype) * size_t(count);
      memcpy(info->rbuf, buf, sz);

      // 受信情報はirecvで発行されているはずのため、waitで削除する
      break;
    }
  }

  return MPI_SUCCESS;
}

/// Blocking receive for a message 
static int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status)
{
  // 逐次なので自身(ランク0)のみ
  int myrank;
  MPI_Comm_rank(comm, &myrank);
  if( source != myrank )
  {
    return MPI_SUCCESS;
  }

  // 同期通信のため、Recvが呼ばれたタイミングで送信情報が存在するはず

  // バッファを検索
  std::vector<CPM_STUBCOMMBUF_INFO*>::iterator it = cpm_StubCommBufInfo.begin();
  for( ;it!=cpm_StubCommBufInfo.end();it++ )
  {
    // tagが同じ送信情報を取得:
    CPM_STUBCOMMBUF_INFO *info = *it;
    if( info->tag == tag && info->sbuf )
    {
      // コピー
      size_t sz = cpm_StubGetDatatypeSize(datatype) * size_t(count);
      memcpy(buf, info->sbuf, sz);

      // 送信情報はisendで発行されているはずのため、waitで削除する
      break;
    }
  }

  return MPI_SUCCESS;
}

/// Begins a nonblocking send 
static int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)
{
  // 逐次なので自身(ランク0)のみ
  int myrank;
  MPI_Comm_rank(comm, &myrank);
  if( dest != myrank )
  {
    return MPI_SUCCESS;
  }

  // request
  *request = cpm_StubGetRequest();

  // 受信情報を格納
  CPM_STUBCOMMBUF_INFO *info = new CPM_STUBCOMMBUF_INFO();
  info->tag = tag;
  info->req = *request;
  info->size = cpm_StubGetDatatypeSize(datatype);
  info->nmemb = size_t(count);
  size_t sz = info->size * info->nmemb;
  info->sbuf = (void*)buf;
  cpm_StubCommBufInfo.push_back(info);

  return MPI_SUCCESS;
}

/// Begins a nonblocking receive 
static int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request *request)
{
  // 逐次なので自身(ランク0)のみ
  int myrank;
  MPI_Comm_rank(comm, &myrank);
  if( source != myrank )
  {
    return MPI_SUCCESS;
  }

  // request
  *request = cpm_StubGetRequest();

  // 送信情報を格納
  CPM_STUBCOMMBUF_INFO *info = new CPM_STUBCOMMBUF_INFO();
  info->tag = tag;
  info->req = *request;
  info->size = cpm_StubGetDatatypeSize(datatype);
  info->nmemb = size_t(count);
  size_t sz = info->size * info->nmemb;
  info->rbuf = (void*)buf;
  cpm_StubCommBufInfo.push_back(info);

  return MPI_SUCCESS;
}

/// Waits for an MPI request to complete 
static int MPI_Wait(MPI_Request *request, MPI_Status *status)
{
  // 逐次なので、この時点でrequestの相手も情報が存在するはず
  // 相手が存在したときにコピーを行う

  // リクエストの情報を検索
  std::vector<CPM_STUBCOMMBUF_INFO*>::iterator it1 = cpm_StubCommBufInfo.begin();
  bool exist = false;
  for( ;it1!=cpm_StubCommBufInfo.end();it1++ )
  {
    if( (*it1)->req = *request )
    {
      exist = true;
      break;
    }
  }
  if( !exist )
  {
    // 存在しなかったのでリターン
    return MPI_SUCCESS;
  }

  // 情報のポインタ
  CPM_STUBCOMMBUF_INFO *infoS = NULL; //受信側
  CPM_STUBCOMMBUF_INFO *infoR = NULL; //送信側
  if( (*it1)->sbuf )
  {
    infoS = *it1;
  }
  else
  {
    infoR = *it1;
  }
  // 自身の情報はバッファから削除
  cpm_StubCommBufInfo.erase(it1);

  // 相手を検索
  std::vector<CPM_STUBCOMMBUF_INFO*>::iterator it2 = cpm_StubCommBufInfo.begin();
  exist = false;
  for( ;it2!=cpm_StubCommBufInfo.end();it2++ )
  {
    if( infoS && (*it2)->rbuf )
    {
      if( infoS->tag == (*it2)->tag )
      {
        infoR = *it2;
        exist = true;
        break;
      }
    }
    if( infoR && (*it2)->sbuf )
    {
      if( (*it1)->tag == (*it2)->tag )
      {
        infoS = *it2;
        exist = true;
        break;
      }
    }
  }

  // 相手の情報もバッファから削除
  if( exist )
  {
    cpm_StubCommBufInfo.erase(it2);
  }

  // 両方そろった場合、コピー
  if( infoS && infoR )
  {
    // コピー
    size_t sz = infoR->size * infoR->nmemb;
    memcpy(infoR->rbuf, infoS->sbuf, sz);
  }

  // 情報は削除
  delete infoS;
  delete infoR;

  return MPI_SUCCESS;
}

/// Waits for all given MPI Requests to complete 
static int MPI_Waitall(int count, MPI_Request array_of_requests[], 
                       MPI_Status array_of_statuses[])
{
  for( int i=0;i<count;i++ )
  {
    MPI_Wait(&array_of_requests[i], &array_of_statuses[i]);
  }
  return MPI_SUCCESS;
}

/// Combines values from all processes and distributes the result back to all processes 
static int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  size_t sz = cpm_StubGetDatatypeSize(datatype);
  if( sz==0 ) return MPI_SUCCESS;
  memcpy(recvbuf, sendbuf, sz*count);
  return MPI_SUCCESS;
}

/// Gathers together values from a group of processes 
static int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm comm)
{
  size_t sz = cpm_StubGetDatatypeSize(recvtype);
  if( sz==0 ) return MPI_SUCCESS;
  memcpy(recvbuf, sendbuf, sz*recvcount);
  return MPI_SUCCESS;
}

/// Gathers data from all tasks and distribute the combined data to all tasks 
static int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm)
{
  size_t sz = cpm_StubGetDatatypeSize(recvtype);
  if( sz==0 ) return MPI_SUCCESS;
  memcpy(recvbuf, sendbuf, sz*recvcount);
  return MPI_SUCCESS;
}

/// Gathers into specified locations from all processes in a group 
static int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, const int *recvcounts, const int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  size_t sz = cpm_StubGetDatatypeSize(recvtype);
  if( sz==0 ) return MPI_SUCCESS;
  char *recvptr = (char*)recvbuf + sz*displs[0];
  memcpy(recvptr, sendbuf, sz*recvcounts[0]);
  return MPI_SUCCESS;
}

/// Gathers data from all tasks and deliver the combined data to all tasks 
static int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, const int *recvcounts, const int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm)
{
  size_t sz = cpm_StubGetDatatypeSize(recvtype);
  if( sz==0 ) return MPI_SUCCESS;
  char *recvptr = (char*)recvbuf + sz*displs[0];
  memcpy(recvptr, sendbuf, sz*recvcounts[0]);
  return MPI_SUCCESS;
}

} // extern "C"

#ifdef __cplusplus
/// namespace for C++ bindings
namespace MPI
{
  /// MPI data types
  enum Datatype
  {
    DATATYPE_NULL      =  0 ///< null
  , CHAR               =  1 ///< char
  , UNSIGNED_CHAR      =  2 ///< unsigned char
  , BYTE               =  3 ///< byte(not support)
  , SHORT              =  4 ///< short
  , UNSIGNED_SHORT     =  5 ///< unsigned short
  , INT                =  6 ///< int
  , UNSIGNED           =  7 ///< unsigned
  , LONG               =  8 ///< long
  , UNSIGNED_LONG      =  9 ///< unsigned long
  , FLOAT              = 10 ///< float
  , DOUBLE             = 11 ///< double
  , LONG_DOUBLE        = 12 ///< long double
  , LONG_LONG_INT      = 13 ///< long long int
  , LONG_LONG          = 13 ///< long long
  , UNSIGNED_LONG_LONG = 51 ///< unsigned long long
  , REAL               = 52 ///< REAL_TYPE
  };

  /// 
  class Intracomm
  {
    public:
    /// Broadcasts a message from the process with rank "root" to all other processes of the communicator 
    int Bcast(void *buffer, int count, Datatype datatype, int root)
    {
      return MPI_Bcast(buffer, count, (MPI_Datatype)datatype, root, 0);
    }
  };

  static Intracomm COMM_WORLD;			///< world communicator
  const int PROC_NULL	= MPI_PROC_NULL;	///< null process

};
#endif /*__cplusplus*/

#endif /*DISABLE_MPI*/
