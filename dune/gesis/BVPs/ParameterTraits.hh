#ifndef DUNE_GESIS_PARAMETER_TRAITS_HH
#define DUNE_GESIS_PARAMETER_TRAITS_HH

namespace Dune {
  namespace Gesis {

    //! Traits class for transport parameters
    template<typename GV,
             typename RF,
             typename IDT,
             typename SDT,
             typename YFG>
    struct GroundwaterParameterTraits {

      typedef GV GridViewType;
      typedef RF RangeFieldType;
      typedef IDT InputDataType;
      typedef SDT SetupDataType;
      typedef YFG YFieldGeneratorType;

      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      typedef typename GV::Grid::ctype DomainFieldType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType, dimDomain - 1 > IntersectionDomainType;

      //! \brief range type
      typedef Dune::FieldVector<RF, GV::dimension> RangeType;

      //! \brief Diffusion tensor type
      typedef Dune::FieldMatrix< RF, GV::dimension, GV::dimension > DiffusionTensorType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

    };






    //! Traits class for transport parameters
    template<typename GV,typename RF,typename DGF,typename IDT,typename SDT>
    struct TransportParameterTraits {
      //! \brief the grid view
      typedef GV GridViewType;

      typedef IDT InputDataType;
      typedef SDT SetupDataType;
      typedef DGF DiscreteGridFunction;

      //! \brief Enum for domain dimension

      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      typedef typename GV::Grid::ctype DomainFieldType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType, dimDomain - 1 > IntersectionDomainType;

      //! \brief Export type for range field
      typedef RF RangeFieldType;

      //! \brief range type
      typedef Dune::FieldVector<RF,GV::dimension> RangeType;

      //! \brief Diffusion tensor type
      typedef Dune::FieldMatrix< RF, GV::dimension, GV::dimension > DispersionTensorType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;
    };


  }  // Gesis
}   // Dune

#endif // DUNE_GESIS_PARAMETER_TRAITS_HH
