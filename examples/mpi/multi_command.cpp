#include <boost/mpi.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <boost/math/tools/promotion.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <iostream>

#include <stan/math.hpp>

struct command {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {}
  virtual void run() const = 0;
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT( command );

struct null_command : public command {
  void run() const {}
};



struct greeter : public command {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(command);
  }
  void run() const {
    boost::mpi::communicator world;
    std::cout << "Hello world from " << world.rank() << "!" << std::endl;
  }
};

BOOST_CLASS_EXPORT(greeter);


struct hard_work {
  double operator()(double a, double b) {
    return a+b;
  }
};

template <typename F>
struct mapped_functor : public command {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(command);
  }
  void run() const {
    boost::mpi::communicator world;
    F fun;
    std::cout << "Job " << world.rank() << " says " << fun(1, world.rank()) << std::endl;
  }
};

BOOST_CLASS_EXPORT(mapped_functor<hard_work>);


struct killer : public command {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    //ar & boost::serialization::base_object<command>(*this);
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(command);
  }
  void run() const {
    boost::mpi::communicator world;
    std::cout << "Terminating worker " << world.rank() << std::endl;
    MPI_Finalize();
    std::exit(0);
  }
};

BOOST_CLASS_EXPORT(killer);


int main(int argc, char* argv[])
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;

  std::cout << "Started worker " << world.rank() << std::endl;
  
  if(world.rank() != 0) {
    while(1) {
      boost::shared_ptr<command> work;

      boost::mpi::broadcast(world, work, 0);

      work->run();
    }
  }

  boost::shared_ptr<command> hello(new greeter);
  
  boost::mpi::broadcast(world, hello, 0);

  hello->run();

  typedef mapped_functor<hard_work> hard_work_t;
  
  boost::shared_ptr<command> work(new hard_work_t);

  boost::mpi::broadcast(world, work, 0);

  boost::shared_ptr<command> quit(new killer);

  boost::mpi::broadcast(world, quit, 0);

  std::cout << "Finally finishing the main worker " << world.rank() << std::endl;
}
