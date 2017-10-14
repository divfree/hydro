#pragma once

#include "std_ref.hpp"
#include "common.hpp"

// ��������� (�����������) � �������� �������
template <class TKey, class TData>
class BinSearchSet
/*
����������: O(log n), ����������: O(n), ��������: O(n)
�������� ������ ���� TData � ������ ���� TKey
��� TKey (����) ������ ���� ����������:
  ����������� �����������, ����������� �� ���������,
  �������� ������������, ��������� ��������� ==, <,>,<=,>=,!=
��� TData (������) ������ ���� ����������:
  ����������� �����������, ����������� �� ���������,
  �������� ������������
*/
{
  TData* data;
  TKey* keys;
public:
  int length;
  BinSearchSet();
  ~BinSearchSet();
  TData* get(TKey); // ���������� ��������� �� ������ � ������ ������
  TData* set(TKey, TData); // ��������� ����� ������ � ������ skey ���
                           // �������� ������������ � ���������� ��������� �� ��
  bool del(TKey); // ������� ������ � ������ ������
  TData value(TKey); // ���������� �������� ������ � ������ ������
  TKey get_key(int i) // ���������� ���� ������ � ������� i
  {
    if (i<length && i>=0) return keys[i]; else return TKey();
  }
  TData get_data(int i) // ���������� �������� ������ � ������� i
  {
    if (i<length && i>=0) return data[i]; else return TData();
  }
  TData* operator ()(TKey key)
  {
    return get(key);
  }
  TData& operator [](TKey key)
  {
    TData* ptr=get(key);
    if(ptr) return *ptr;
    else throw string("'"+key+"' undefined");
  }
  bool exist(TKey key)
  {
    return (get(key)!=0);
  }
};

template <class TKey, class TData>
TData* BinSearchSet<TKey, TData>::set(TKey skey, TData sdata)
{
  int L=0, R=length-1, M=0;
  while (L<=R)
  {
    M=(L+R)/2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1;
        else
        {
          data[M]=sdata;
          return &data[M];
        }
  }
  int i=max(0, M-1);
  while (i<length && keys[i]<skey) i++;
  TData* tdata = new TData[length+1];
  TKey* tkeys = new TKey[length+1];
  for (int j=0; j<i; j++)
  {
    tdata[j]=data[j];
    tkeys[j]=keys[j];
  }
  tdata[i]=sdata;
  tkeys[i]=skey;
  for (int j=i+1; j<length+1; j++)
  {
    tdata[j]=data[j-1];
    tkeys[j]=keys[j-1];
  }
  length++;
  delete [] data;
  delete [] keys;
  data = tdata;
  keys = tkeys;
  return &data[i];
}
template <class TKey, class TData>
BinSearchSet<TKey, TData>::BinSearchSet()
{
  length = 0;
  data = 0;
  keys = 0;
}
template <class TKey, class TData>
BinSearchSet<TKey, TData>::~BinSearchSet()
{
  delete [] data;
  delete [] keys;
}
template <class TKey, class TData>
TData* BinSearchSet<TKey, TData>::get(TKey skey)
{
  int L=0, R=length-1;
  while (L<=R)
  {
    int M = (L+R) / 2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1; else
        return &data[M];
  }
  return 0;
}
template <class TKey, class TData>
TData BinSearchSet<TKey, TData>::value(TKey skey)
{
  int L=0, R=length-1;
  while (L<=R)
  {
    int M = (L+R) / 2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1; else
        return data[M];
  }
  return TData();
}
template <class TKey, class TData>
bool BinSearchSet<TKey, TData>::del(TKey skey)
{
  int L=0, R=length-1;
  while (L<=R)
  {
    int M = (L+R) / 2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1;
        else
          if (length==1)
          {
            delete [] data;
            delete [] keys;
            keys=0;
            data=0;
            length=0;
            return true;
          } else
          {
            TData* tdata=new TData[length-1];
            TKey* tkeys=new TKey[length-1];
            for (int j=0; j<M; j++)
            {
              tdata[j]=data[j];
              tkeys[j]=keys[j];
            }
            for (int j=M; j<length-1; j++)
            {
              tdata[j]=data[j+1];
              tkeys[j]=keys[j+1];
            }
            length--;
            delete [] data;
            delete [] keys;
            data=tdata;
            keys=tkeys;
            return true;
          }
  }
  return false;
}

// ����������� �����
// �������� -- ��������� (�����������) � ������
template <class TKey>
class TParameter
{
public:
  BinSearchSet<string, int> PK; // Possible Keys
  void import_PK(int length, TKey* list)
  {
    for(int i=0; i<length; i++)
      PK.set(list[i],0);
  }
  virtual ~TParameter() {}
  virtual int get_length() const = 0;
  virtual bool exist(TKey key) = 0;
  virtual bool del(TKey) = 0; // ������� ������ � ������ ������
  virtual TKey get_key(int i) const = 0; // ���������� ���� ������ � ������� i
  virtual void read_data_by_key(istream& in, TKey key) = 0;
  virtual void write_data_by_index(ostream& out, int i) = 0;
  virtual bool write_data_by_key(ostream& out, TKey key) = 0;
};



// ��������� (�����������) � �������� ������� � ������������ �����/������
template <class TKey, class TData>
class BinSearchSet2 : public TParameter<TKey>
/*
����������: O(log n), ����������: O(n), ��������: O(n)
�������� ������ ���� TData � ������ ���� TKey
��� TKey (����) ������ ���� ����������:
  ����������� �����������, ����������� �� ���������,
  �������� ������������, ��������� ��������� ==, <,>,<=,>=,!=
��� TData (������) ������ ���� ����������:
  ����������� �����������, ����������� �� ���������,
  �������� ������������
  operator<<(ostream)
  operator>>(istream)
*/
{
  TData* data;
  TKey* keys;
public:
  int length;
  BinSearchSet2();
  ~BinSearchSet2();
  TData* get(TKey); // ���������� ��������� �� ������ � ������ ������
  TData* set(TKey, const TData&); // ��������� ����� ������ � ������ skey ���
                           // �������� ������������ � ���������� ��������� �� ��
  bool del(TKey); // ������� ������ � ������ ������
  TData value(TKey); // ���������� �������� ������ � ������ ������
  TKey get_key(int i) const // ���������� ���� ������ � ������� i
  {
    if (i<length && i>=0) return keys[i]; else return TKey();
  }
  TData get_data(int i) const // ���������� �������� ������ � ������� i
  {
    if (i<length && i>=0) return data[i]; else return TData();
  }
  TData* operator ()(TKey key)
  {
    return get(key);
  }
  TData& operator [](TKey key)
  {
    TData* ptr=get(key);
    if(ptr) return *ptr;
    else throw string("'"+key+"' undefined");
  }
  int get_length() const
  {
    return length;
  }
  bool exist(TKey key)
  {
    return (get(key)!=0);
  }
  void read_data_by_key(istream& in, TKey key)
  {
    TData t;
    read_data(in, t);
    set(key, t);
  }
  void write_data_by_index(ostream& out, int i)
  {
    out<<data[i];
  }
  bool write_data_by_key(ostream& out, TKey key)
  {
    TData* ptr=get(key);
    if(ptr) out<<*ptr;
    return (ptr!=0);
  }
};

template <class TKey, class TData>
TData* BinSearchSet2<TKey, TData>::set(TKey skey, const TData& sdata)
{
  int L=0, R=length-1, M=0;
  while (L<=R)
  {
    M=(L+R)/2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1;
        else
        {
          data[M]=sdata;
          return &data[M];
        }
  }
  int i=max(0, M-1);
  while (i<length && keys[i]<skey) i++;
  TData* tdata = new TData[length+1];
  TKey* tkeys = new TKey[length+1];
  for (int j=0; j<i; j++)
  {
    tdata[j]=data[j];
    tkeys[j]=keys[j];
  }
  tdata[i]=sdata;
  tkeys[i]=skey;
  for (int j=i+1; j<length+1; j++)
  {
    tdata[j]=data[j-1];
    tkeys[j]=keys[j-1];
  }
  length++;
  delete [] data;
  delete [] keys;
  data = tdata;
  keys = tkeys;
  return &data[i];
}
template <class TKey, class TData>
BinSearchSet2<TKey, TData>::BinSearchSet2()
{
  length = 0;
  data = 0;
  keys = 0;
}
template <class TKey, class TData>
BinSearchSet2<TKey, TData>::~BinSearchSet2()
{
  delete [] data;
  delete [] keys;
}
template <class TKey, class TData>
TData* BinSearchSet2<TKey, TData>::get(TKey skey)
{
  int L=0, R=length-1;
  while (L<=R)
  {
    int M = (L+R) / 2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1; else
        return &data[M];
  }
  return 0;
}
template <class TKey, class TData>
TData BinSearchSet2<TKey, TData>::value(TKey skey)
{
  int L=0, R=length-1;
  while (L<=R)
  {
    int M = (L+R) / 2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1; else
        return data[M];
  }
  return TData();
}
template <class TKey, class TData>
bool BinSearchSet2<TKey, TData>::del(TKey skey)
{
  int L=0, R=length-1;
  while (L<=R)
  {
    int M = (L+R) / 2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1;
        else
          if (length==1)
          {
            delete [] data;
            delete [] keys;
            keys=0;
            data=0;
            length=0;
            return true;
          } else
          {
            TData* tdata=new TData[length-1];
            TKey* tkeys=new TKey[length-1];
            for (int j=0; j<M; j++)
            {
              tdata[j]=data[j];
              tkeys[j]=keys[j];
            }
            for (int j=M; j<length-1; j++)
            {
              tdata[j]=data[j+1];
              tkeys[j]=keys[j+1];
            }
            length--;
            delete [] data;
            delete [] keys;
            data=tdata;
            keys=tkeys;
            return true;
          }
  }
  return false;
}


// ��������� � �������� �������, �������� ��������� �� �������
template <class TKey, class TData>
class BinSearchSet_ptr
/*
����������: O(log n), ����������: O(n), ��������: O(n)
�������� ������ ���� TData � ������ ���� TKey
��� TKey (����) ������ ���� ����������:
  ����������� �����������, ����������� �� ���������,
  �������� ������������, ��������� ��������� ==, <,>,<=,>=,!=
��� TData (������) ������ ���� ����������:
  ����������� �� ���������
*/
{
  TData** data;
  TKey* keys;
public:
  int length;
  BinSearchSet_ptr();
  ~BinSearchSet_ptr();
  TData* get(TKey); // ���������� ��������� �� ������ � ������ ������
  TData* set(TKey); // ��������� ����� ������
  bool del(TKey); // ������� ������
  TKey get_key(int i) // ���������� ���� ������ � ������� i
  {
    if (i<length && i>=0) return keys[i]; else return TKey();
  }
  TData* get_ptr(int i) // ���������� ��������� �� �������� ������ � ������� i
  {
    if (i<length && i>=0) return data[i]; else return 0;
  }
  TData* operator ()(TKey key)
  {
    return get(key);
  }
  TData& operator [](TKey key)
  {
    return *get(key);
  }
};

template <class TKey, class TData>
TData* BinSearchSet_ptr<TKey, TData>::set(TKey skey)
{
  int L=0, R=length-1, M=0;
  while (L<=R)
  {
    M=(L+R)/2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1;
        else
        {
          // the entry exists
          return data[M];
        }
  }
  // new entry
  int i=max(0, M-1);
  while (i<length && keys[i]<skey) i++;
  TData** tdata = new TData*[length+1];
  TKey* tkeys = new TKey[length+1];
  for (int j=0; j<i; j++)
  {
    tdata[j]=data[j];
    tkeys[j]=keys[j];
  }
  tdata[i]=new TData;
  tkeys[i]=skey;
  for (int j=i+1; j<length+1; j++)
  {
    tdata[j]=data[j-1];
    tkeys[j]=keys[j-1];
  }
  length++;
  delete [] data;
  delete [] keys;
  data = tdata;
  keys = tkeys;
  return data[i];
}
template <class TKey, class TData>
BinSearchSet_ptr<TKey, TData>::BinSearchSet_ptr()
{
  length = 0;
  data = 0;
  keys = 0;
}
template <class TKey, class TData>
BinSearchSet_ptr<TKey, TData>::~BinSearchSet_ptr()
{
  for(int i=0; i<length; i++) delete data[i];
  delete [] data;
  delete [] keys;
}
template <class TKey, class TData>
TData* BinSearchSet_ptr<TKey, TData>::get(TKey skey)
{
  int L=0, R=length-1;
  while (L<=R)
  {
    int M = (L+R) / 2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1; else
        return data[M];
  }
  return 0;
}
template <class TKey, class TData>
bool BinSearchSet_ptr<TKey, TData>::del(TKey skey)
{
  int L=0, R=length-1;
  while (L<=R)
  {
    int M = (L+R) / 2;
    if (keys[M]>skey) R=M-1; else
      if (keys[M]<skey) L=M+1;
        else
        {
          delete data[M];
          if (length==1)
          {
            delete [] data;
            delete [] keys;
            keys=0;
            data=0;
            length=0;
            return true;
          } else
          {
            TData** tdata=new TData*[length-1];
            TKey* tkeys=new TKey[length-1];
            for (int j=0; j<M; j++)
            {
              tdata[j]=data[j];
              tkeys[j]=keys[j];
            }
            for (int j=M; j<length-1; j++)
            {
              tdata[j]=data[j+1];
              tkeys[j]=keys[j+1];
            }
            length--;
            delete [] data;
            delete [] keys;
            data=tdata;
            keys=tkeys;
            return true;
          }
      }
  }
  return false;
}

