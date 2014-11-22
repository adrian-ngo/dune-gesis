#ifndef DUNE_GEOINVERSION_REDISTRIBUTE_DATAHANDLE_HH
#define DUNE_GEOINVERSION_REDISTRIBUTE_DATAHANDLE_HH

//! Data handle for load_balance()
template<typename GridType,typename VCType>
class RedistributeDataHandle :
  public Dune::CommDataHandleIF<RedistributeDataHandle<GridType,VCType>, int>
{

  typedef typename GridType::LocalIdSet IdSet;

  typedef typename IdSet::IdType Id;
  const GridType &grid;
  std::map<Id,REAL> &elementTags;

public:
  RedistributeDataHandle(const GridType &grid_,
                         std::map<Id,REAL> &elementTags_
                         ) :
    grid(grid_),
    elementTags(elementTags_)
  { };

  bool contains(int dim, int codim) const
  { return codim == 0; }

  bool fixedsize (int dim, int codim) const {
    return false;
  }

  template<class Entity>
  std::size_t size (const Entity &e) const {
    DUNE_THROW(Dune::RangeError, "Nothing needs to be communicated for codim "
               << Entity::codimension << ", so this method (std::size_t " <<
               Dune::className(*this) << "::size(const " <<
               Dune::className<Entity>() << "&) const) should never be "
               "called.");
  }

  std::size_t size (const typename GridType::template Codim<0>::Entity &e) const
  {
    return elementTags.count(grid.localIdSet().id(e));
  }


  template<class MessageBuffer, class Entity>
  void gather(MessageBuffer &buff, const Entity &e) const {
    DUNE_THROW(Dune::RangeError, "Nothing needs to be communicated for codim "
               << Entity::codimension << ", so this method (void "
               << Dune::className(*this) << "::gather(" <<
               Dune::className<MessageBuffer>() << "&, const " <<
               Dune::className<Entity>() << "&) const) should never be "
               "called.");
  }

  /* Sender */
  template<class MessageBuffer>
  void gather(MessageBuffer &buff,
              const typename GridType::template Codim<0>::Entity &e) const
  {

    typename std::map<Id,REAL>::const_iterator it =
      elementTags.find(grid.localIdSet().id(e));
    if(it != elementTags.end())
      buff.write(it->second);
  }

  template<class MessageBuffer, class Entity>
  void scatter(MessageBuffer &buff, const Entity &e, std::size_t n) {
    DUNE_THROW(Dune::RangeError, "Nothing needs to be communicated for codim "
               << Entity::codimension << ", so this method (void " <<
               Dune::className(*this) << "::scatter(" <<
               Dune::className<MessageBuffer>() << "&, const " <<
               Dune::className<Entity>() << "&, std::size_t)) should never be "
               "called.");
  }

  /* Receiver */
  template<class MessageBuffer>
  void scatter(MessageBuffer &buff,
               const typename GridType::template Codim<0>::Entity &e,
               std::size_t n)
  {
    switch(n) {
    case 0: break;
    case 1:
      buff.read(elementTags[grid.localIdSet().id(e)]);
      break;
    default:
      std::cout << std::endl;
      std::cout << "WARNING: grid.localIdSet().id(e) = " << grid.localIdSet().id(e) << std::endl;
      std::cout << "WARNING: grid.globalIdSet().id(e) = " << grid.globalIdSet().id(e) << std::endl;
      std::cout << "WARNING: scatter n=" << n << std::endl;
      DUNE_THROW(Dune::RangeError, "At most one data item may be "
                 "communicated!");
    }
  }

  void compress() { }
};



#endif // DUNE_GEOINVERSION_REDISTRIBUTE_DATAHANDLE_HH
