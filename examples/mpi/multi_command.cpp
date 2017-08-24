#include <iostream>
#include <tuple>
#include <string>

#include <boost/mpi.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <boost/math/tools/promotion.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <stan/math.hpp>

// assume the user has defined something like this in data or
// transformed data (in that order!)
// real foo:
// real bar[3];
// int ind[2];
typedef std::tuple<double, std::vector<double>, std::vector<int> > static_data;

struct command {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {}
  virtual void run(const static_data& data) const = 0;
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT( command );

struct null_command : public command {
  void run(const static_data& data) const {}
};



struct greeter : public command {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(command);
  }
  void run(const static_data& data) const {
    boost::mpi::communicator world;
    std::cout << "Hello world from " << world.rank() << "!" << std::endl;
  }
};

BOOST_CLASS_EXPORT(greeter);

// this struct glues together what user function to call with what
// element from the static_data tuple
template <typename F, int A1, int A2>
struct binary_data_function {
  double operator()(double param, const static_data& data) {
    return(F()(param, std::get<A1>(data), std::get<A2>(data)));
  }
};

// the actual work to be done on the workers; So the user just defines
// as usual a normal function which is turned into a usual functor
// object.
struct hard_work {
  double operator()(double param, const std::vector<double>& bar, const std::vector<int>& ind) {
    return param + bar[0];
  }
};

template <typename F>
struct mapped_functor : public command {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(command);
  }
  void run(const static_data& data) const {
    boost::mpi::communicator world;
    // recieve some parameters here... fix it for now
    double param = world.rank();
    F fun;
    std::cout << "Job " << world.rank() << " says " << fun(param, data) << std::endl;
  }
};

typedef binary_data_function<hard_work,1,2> user_data_function;
BOOST_CLASS_EXPORT(mapped_functor<user_data_function>);


struct killer : public command {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    //ar & boost::serialization::base_object<command>(*this);
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(command);
  }
  void run(const static_data& data) const {
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

  double foo = 10.;
  std::vector<double> bar(3, 5.);
  std::vector<int> ind(2, 1);
  static_data user_data = make_tuple(foo, bar, ind);

  std::cout << "Started worker " << world.rank() << std::endl;
  
  if(world.rank() != 0) {
    while(1) {
      boost::shared_ptr<command> work;

      boost::mpi::broadcast(world, work, 0);

      work->run(user_data);
    }
  }

  boost::shared_ptr<command> hello(new greeter);
  
  boost::mpi::broadcast(world, hello, 0);

  hello->run(user_data);

  typedef mapped_functor<user_data_function> user_work_t;
  
  boost::shared_ptr<command> work(new user_work_t);

  boost::mpi::broadcast(world, work, 0);

  boost::shared_ptr<command> quit(new killer);

  boost::mpi::broadcast(world, quit, 0);

  std::cout << "Finally finishing the main worker " << world.rank() << std::endl;
}
