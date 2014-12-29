/* 
 * File:   FieldData.hh
 * Author: A. Ngo (09/2014)
 *
 * Data structure used to read the input data for the random field to be generated.
 * Corresponds to "yfield.cfg".
 *
 */

#ifndef FIELD_DATA_HH
#define	FIELD_DATA_HH

namespace Dune {
  namespace Gesis {

    class FieldData {
    public:
      int dim;
      Vector<REAL> extensions;
      Vector<unsigned int> nCells;
      Vector<unsigned int> nCellsExt;
      Vector<REAL> gridsizes;
      Vector<REAL> correlations;
      REAL beta;
      REAL fieldVariance;
      std::string variogramModel;
      REAL embeddingFactor;
      bool newField;
      bool newEV;
      bool showEV;
      int seed;
      std::string location;

      FieldData( const int dim_ ){
        dim = dim_;
        extensions.resize(dim);
        nCells.resize(dim);
        nCellsExt.resize(dim);
        gridsizes.resize(dim);
        correlations.resize(dim);
      };

      template<typename IDT> 
      FieldData( const IDT& inputdata_ ){
        dim = inputdata_.domain_data.dim;
        extensions.resize(dim);
        nCells.resize(dim);
        nCellsExt.resize(dim);
        gridsizes.resize(dim);
        correlations.resize(dim);
        retrieve( inputdata_ );
        extend_domain();
      };


      template<typename IDT> 
      void retrieve(const IDT& inputdata){

        extensions = inputdata.domain_data.extensions;
        nCells     = inputdata.domain_data.nCells;

        for(int i=0;i<dim;++i)
          gridsizes[i] = extensions[i] / REAL( nCells[i] );

        correlations = inputdata.yfield_properties.zones[0].correlation_lambda;

        beta            = inputdata.yfield_properties.zones[0].beta;
        fieldVariance   = inputdata.yfield_properties.zones[0].variance;
        embeddingFactor = inputdata.yfield_properties.zones[0].embedding_factor;
        variogramModel  = inputdata.yfield_properties.zones[0].model;
        newField        = inputdata.problem_types.new_YField;
        newEV           = inputdata.problem_types.new_Eigenvalues;
        showEV          = (inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL) ? true : false;
        seed            = inputdata.yfield_properties.random_seed;
        
      };


      void extend_domain() {
        // compute the size of the extended domain required for circulant embedding
        for (int i = 0; i < dim; i++) {
          REAL b = extensions[i];
          REAL a = 2.0 * b;
          REAL f = embeddingFactor;
          REAL c = correlations[i];
          REAL d = gridsizes[i];
          nCellsExt[i] = std::ceil( std::max( a , b + f * c ) / d );
          //nCells_ExtendedDomain[i] = 2 * ( fielddata.nCells[i] - 1 );
        }
      };



      template<typename CFG> 
      void read(const CFG& configuration){

        extensions[0] = configuration.template get<REAL>("domain.Lx");
        extensions[1] = configuration.template get<REAL>("domain.Ly");
        if( dim == 3 )
          extensions[2] = configuration.template get<REAL>("domain.Lz");

        nCells[0] = configuration.template get<int>("grid.nx");
        nCells[1] = configuration.template get<int>("grid.ny");
        if( dim == 3 )
          nCells[2] = configuration.template get<int>("grid.nz");

        for(int i=0;i<dim;++i)
          gridsizes[i] = extensions[i] / REAL( nCells[i] );

        correlations[0] = configuration.template get<REAL>("field.lx");
        correlations[1] = configuration.template get<REAL>("field.ly");
        if( dim == 3 )
          correlations[2] = configuration.template get<REAL>("field.lz");

        beta = configuration.template get<REAL>("field.beta");
        fieldVariance = configuration.template get<REAL>("field.variance");
        embeddingFactor = configuration.template get<REAL>("field.EF");
        variogramModel = configuration.template get<std::string>("field.VM");
        newField = ("true"==configuration.template get<std::string>("field.newField"));
        newEV = ("true"==configuration.template get<std::string>("field.newEV"));
        showEV = ("true"==configuration.template get<std::string>("field.showEV"));
        seed = configuration.template get<int>("field.seed");
        location = configuration.template get<std::string>("field.location");
        
      };

      void printInfos(){
        std::cout << "Dimension           : " << dim << std::endl;
        std::cout << "Domain              : " << General::Vector2CSV( extensions ) << std::endl;
        std::cout << "Grid                : " << General::Vector2CSV( nCells ) << std::endl;
        std::cout << "cell-size           : " << General::Vector2CSV( gridsizes ) << std::endl;
        std::cout << "Variogram model     : " << variogramModel << std::endl;
        std::cout << "Correlation-lengths : " << General::Vector2CSV( correlations ) << std::endl;
        std::cout << "Field trend         : " << beta << std::endl;
        std::cout << "Field variance      : " << fieldVariance << std::endl;
        std::cout << "embedding factor    : " << embeddingFactor << std::endl;
        std::cout << "seed                : " << seed << std::endl;
        std::cout << "Extended Grid       : " << General::Vector2CSV( nCellsExt ) << std::endl;
      };
    };

  }
}

#endif	/* FIELD_DATA_HH */
