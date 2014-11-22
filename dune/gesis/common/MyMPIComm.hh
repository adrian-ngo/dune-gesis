// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 * File:   MyMPIComm.hh
 * Author: ronnie
 *
 * similar to the MPI helper. needed to get several MPI communicators.
 * --> used for the parallel calculation of sensitivities!
 *
 *
 */
#ifndef DUNE_GESIS_MYMPICOMM_HH
#define DUNE_GESIS_MYMPICOMM_HH

// get the logger running
extern CLogfile logger;

/*
 * MyMPIComm class definition
 */
class MyMPIComm {
public:
  //1. constructor (empty)
  MyMPIComm():my_group(false){size=-1;rank=-1;};

  //2. constructor
  MyMPIComm(MPI_Comm _comm, MPI_Group _group,bool _my_group=false):comm(_comm),group(_group),my_group(_my_group){
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_world);
    MPI_Comm_size(MPI_COMM_WORLD,&size_world);
    if(my_group){
      MPI_Comm_rank(comm,&rank);
      MPI_Comm_size(comm,&size);
    }else{
      size=-1;
      rank=-1;
    }
  };
  //destructor ->not needed

  /*
   * functions
   */

  // set the given variables!
  void set(MPI_Comm _comm,MPI_Group _group, bool _my_group=false){
    comm=_comm;
    group=_group;
    my_group=_my_group;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_world);
    MPI_Comm_size(MPI_COMM_WORLD,&size_world);
    if(my_group){
      MPI_Comm_rank(comm,&rank);
      MPI_Comm_size(comm,&size);
    }
  };

  //update the rank and size for MPI_COMM_WORLD and the stored MPI communicator
  // IMPORTANT: my_group needs to be set via "set_I_am_in()" by hand before calling update!!!!
  void update(){
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_world);
    MPI_Comm_size(MPI_COMM_WORLD,&size_world);
    if(my_group){
      MPI_Comm_rank(comm,&rank);
      MPI_Comm_size(comm,&size);
    }
  };

  // set and get routines
  void set_I_am_in(){my_group=true;}
  int get_rank() const {return rank;};
  int get_size() const {return size;};
  int get_rank_world() const {return rank_world;};
  int get_size_world() const {return size_world;};
  bool I_am_in() const {return my_group;};

  // get the stored communicator
  MPI_Comm get_comm() const {return comm;};
  MPI_Comm* get_comm_adress(){return &comm;}

  // get the stored group
  MPI_Group get_group() const {return group;};
  MPI_Group* get_group_adress(){return &group;}

  // print some info: for debugging usefull
  void print_info() const {
    logger<<"MyMPIComm INFO: "<<std::endl;
    logger<<"rank_world = "<<rank_world<<std::endl;
    logger<<"size_world = "<<size_world<<std::endl;
    if(rank>-1){
      logger<<"rank = "<<rank<<std::endl;
      logger<<"size = "<<size<<std::endl;
    }else{
      logger<<"I am global P"<<rank_world<<" and NOT included in this pool! (no rank, size here!)"<<std::endl;
    }
    logger<<"MyMPIComm INFO END! "<<std::endl;

  };

  /*
   * private variables!
   */
private:
  int rank; // rank of the process in the group
  int size; // size of the group
  int rank_world; // rank of the processor in MPI_COMM_WORLD
  int size_world; // size of MPI_COMM_WORLD
  MPI_Comm comm; // communicator of the group
  MPI_Group group; // the communicator group!
  bool my_group;  // flag wether the processor is in the group or not
};

#endif // DUNE_GESIS_MYMPICOMM_HH
