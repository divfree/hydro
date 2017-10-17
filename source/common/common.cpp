#include "common.hpp"

string get_time(string fmt)
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer [80];

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  strftime (buffer,80,fmt.c_str(),timeinfo);
  return buffer;
}

string get_string_with_ws(istream& in)
{
  char c=0;
  string res="";
  in>>c;
  if(c!=0)
  {
    if(c=='"')
    {
      getline(in,res,'"');
    } else
    if(in.peek()==' ')
    {
      res=c;
    } else
    {
      string add="";
      in>>add;
      res=c+add;
    }
  }
  return res;
}

string get_current_directory()
{
  //boost::filesystem::path full_path( boost::filesystem::current_path() );
  //return full_path.string();
  assert(false);
  return "";
}

void set_current_directory(string path)
{
  //boost::filesystem::current_path(path);
  assert(false);
}

string IntToStr(int a)
{
  ostringstream stream;
  stream << a;
  return stream.str();
}
string IntToStr(size_t a)
{
  ostringstream stream;
  stream << a;
  return stream.str();
}
string DoubleToStr(double a, int precision)
{
  ostringstream stream;
  stream.precision(precision);
  stream << a;
  return stream.str();
}

double StrToDouble(string s)
{
  istringstream stream;
  stream.str(s);
  double a;
  stream >> a;
  return a;
}

string point_to_string(int i, int j)
{
  return "("+IntToStr(i)+","+IntToStr(j)+")";
}

bool file_exists(string filename)
{
  ifstream ifile(filename);
  if(ifile.good())
  {
    ifile.close();
    return true;
  }
  return false;
}

bool directory_exists(string dirName)
{
  return true;
  /*
  if (boost::filesystem::exists(dirName.c_str()) != 0)
  {
    if (errno == ENOENT || errno == ENOTDIR)
    {
       return false;
    }
    return false;
  }
  return true;
  */
}

bool ecast(bool* ptr)
{
  return (ptr!=0 && *ptr);
}

string ecast(string* ptr)
{
  if(ptr) return *ptr; else return "";
}

int ecast(int* ptr, int def)
{
  if(ptr) return *ptr; else return def;
}

double ecast(double* ptr, double def)
{
  if(ptr) return *ptr; else return def;
}

double sqr(double a)
{
  return a*a;
}

int simple_pow(int x, int y)
{
  int res=1;
  for(int i=0; i<y; i++)
  {
    res*=x;
  }
  return res;
}

// retrun k : b=10^k-1, b>=a, k->min
int sup_power10_k(int a)
{
  int k=0;
  int b=1;
  while(a>b-1)
  {
    b = b*10;
    k++;
  }
  return k;
}

// retrun k : b=10^k, b>=a, k->min
int sup_power10_k(double a)
{
  int k=0;
  double b=1;
  while(a>b)
  {
    b = b*10;
    k++;
  }
  return k;
}

// retrun b : b=10^k-1, b>=a, k->min
int sup_power10(int a)
{
  int b=1;
  while(a>b-1)
  {
    b = b*10;
  }
  return b;
}

// round a to nearest number from {n*10^(-k), n--integer}, k>=0
double round_to_exponent(double a, int k)
{
  double b=1;
  for(int i=0; i<k; i++) b*=10;

  return double(round(a*b))/b;
}

// round a with specified precision
double round_to_precision(double a, int precision)
{
  stringstream buf;
  buf.precision(precision);
  buf<<a;
  double res;
  buf>>res;
  return res;
}

string IntToStr_lz(int a, int min_digits)
{
  int b=abs(a);
  int real_digits=max(sup_power10_k(b),1);
  int add_digits=max(0, min_digits-real_digits);

  string res;
  if(a<0) res="-"; else res="";
  for(int i=0; i<add_digits; i++) res+="0";
  res+=IntToStr(b);
  return res;
}

string seconds_to_hms(double sec)
{
  double rest=sec;
  int hours=int(sec/3600.0); rest-=hours*3600.0;
  int minutes=int(rest/60.0); rest-=minutes*60.0;
  int seconds=int(rest); rest-=seconds;
  int milliseconds=int(rest*1000);

  string res;
  res=IntToStr_lz(hours,2)+":"+IntToStr_lz(minutes,2)+":"+IntToStr_lz(seconds,2)+"."+IntToStr_lz(milliseconds,3);

  return res;
}

double max(double a1, double a2, double a3)
{
  return max(a1, max(a2,a3));
}

int round05(double a)
{
  return int(a+0.5);
}

/*int round(double a)
{
  return boost::math::iround(a);
}*/

int myround(double a)
{
  return round05(a);
}

bool in_rect(int i, int j, int i1, int j1, int i2, int j2)
{
  return (i>=i1 && i<=i2 && j>=j1 && j<=j2);
}

bool in_rect(int i, int j, const TRect& rect)
{
  return (i>=rect.ai && i<=rect.bi && j>=rect.aj && j<=rect.bj);
}

int sysinfo::virtual_usage_kb()
{ //Note: this value is in KB!
  ifstream file("/proc/self/status");

  string parname="";
  bool found=false;

  while(file.good() && !file.eof())
  {
    file>>parname;
    if(parname=="VmSize:")
    {
      found=true;
      break;
    }
  }

  int a=0;
  if(found)
  {
    file>>a;
  }

  return a;
}

int sysinfo::physical_usage_kb()
{ //Note: this value is in KB!
  ifstream file("/proc/self/status");

  string parname="";
  bool found=false;

  while(file.good() && !file.eof())
  {
    file>>parname;
    if(parname=="VmRSS:")
    {
      found=true;
      break;
    }
  }

  int a=0;
  if(found)
  {
    file>>a;
  }

  return a;
}

int sysinfo::threads()
{ //Note: this value is in KB!
  ifstream file("/proc/self/status");

  string parname="";
  bool found=false;

  while(file.good() && !file.eof())
  {
    file>>parname;
    if(parname=="Threads:")
    {
      found=true;
      break;
    }
  }

  int a=0;
  if(found)
  {
    file>>a;
  }

  return a;
}

void create_directory(string dirname)
{
  /*
  using namespace boost::filesystem;
	path dir(dirname.c_str());
  try
  {
    create_directories(dir);
  }
  catch (const filesystem_error& ex)
  {
    cout << ex.what() << '\n';
      throw string("Can't create the directory: '"+dirname+"'");
  }
  */
  assert(false);
}

void remove_file(string filename)
{
  std::remove(filename.c_str());
}

double pow3(double x)
{
  return x*x*x;
}

double pow4(double x)
{
  return x*x*x*x;
}

double bilinear_interpolation(double x, double y, double x00, double y00, double x11, double y11, double a00, double a10, double a11, double a01)
{
  double b00=a00*(x11-x)*(y11-y);
  double b10=a10*(x-x00)*(y11-y);
  double b11=a11*(x-x00)*(y-y00);
  double b01=a01*(x11-x)*(y-y00);
  return (b00+b10+b11+b01)/((x11-x00)*(y11-y00));
}

double linear_interpolation(double x, double x0, double x1, double a0, double a1)
{
  return ((x-x0)*a1+(x1-x)*a0)/(x1-x0);
}

double distance2(double x1, double y1, double x2, double y2)
{
  return sqrt(sqr(x2-x1)+sqr(y2-y1));
}


void read_data(istream& in, string& data)
{
  data=get_string_with_ws(in);
}
