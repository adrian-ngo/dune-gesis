/*
 * File:   logfile.hh
 * Author: ngo
 *
 * Modified on Sept.10, 2010
 *
 * Usage example inside C++:
 *
 * global version:
 * ===============
 * // in main.cc:
 * CLogfile logger; // declare this globally, outside main()
 *
 * // inside main():
 * std::stringstream slogfilename;
 * slogfilename << "Logfile_P" << helper.rank() << ".log";
 * logger.init( slogfilename.str() );
 * logger << "your output i = " << i << std::endl;
 *
 * // in external header:
 * extern CLogfile logger; // declare this globally, outside functions/classes
 *
 * // inside your functions/classes, you can use it normally:
 * logger << "your output i = " << i << std::endl;
 *
 *
 * // obsolete (local) version:
 * std::stringstream slogfilename;
 * slogfilename << "Logfile_P" << helper.rank() << ".log";
 * CLogfile logger( helper, slogfilename.str() );
 * logger << "your output i = " << i << std::endl;
 *
 */


#ifndef DUNE_GESIS_LOGFILE_HH
#define	DUNE_GESIS_LOGFILE_HH

class CLogfile : public std::ofstream
{
private:
  std::string filename;
public:
  CLogfile(){};

  void init( const std::string filename_ ){
    filename = filename_;
    std::cout << "Open logger file " << filename << std::endl;
    this->open( filename.c_str() );
  };
  
  ~CLogfile(){
    std::cout << "Close logger file " << filename << std::endl;
    this->close();
  };
  
};






#endif // DUNE_GESIS_LOGFILE_HH
