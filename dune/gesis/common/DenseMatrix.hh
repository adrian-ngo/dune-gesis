/*
 * File:   DenseMatrix.hh
 * Author: A. Ngo
 *
 * 2010-2014
 */

#ifndef DUNE_GESIS_DENSEMATRIX_HH
#define	DUNE_GESIS_DENSEMATRIX_HH

#include <vector>
#include <cmath>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "Vector.hh"

#include "dune/gesis/common/io/IO_routines.hh"

namespace Dune {
  namespace Gesis {

    template<typename ComponentType>
    class DenseMatrix
    {
    private:
      Vector<ComponentType> data;
      int rows;
      int cols;
      Vector<ComponentType> row_sums; // for row equilibration
      bool row_equilibrated;

      inline ComponentType & at(const int row, const int col)
      {
        return data[row * cols + col];
      }

      inline const ComponentType & at(const int row, const int col) const
      {
        return data[row * cols + col];
      }

      inline int
      pivot_index(const DenseMatrix & A, const int & row, const int & col, const std::vector<int> & r) const
      {
        //    ComponentType max = A(r[row],col);
        // Korrektur 1: A(r[row],col) darf auch negativ sein. max sollte nur ungleich 0 sein!

        ComponentType max = std::abs(A(r[row],col));

        int imax = row;

        for(int i=row+1; i<rows; ++i)
          {

            //		std::cout << " max = " << max << std::endl;
            //		std::cout << " A(r[row],col) = " << abs(A(r[row],col)) << std::endl;
            //		std::cout << " A(r[i],col)   = " << abs(A(r[i],col)) << std::endl;

            if( std::abs(A(r[i],col)) > std::abs(A(r[row],col) ) )
              {
                //			max = A(r[i],col);
                // Korrektur 1: A(r[row],col) darf auch negativ sein. max sollte nur ungleich 0 sein!
                max = std::abs(A(r[i],col)); //A(r[i],col);
                imax = i;
              }
          }
        if( max < -1E-12 )
          std::cout << "Warning: DenseMatrix::pivot_index(): max < 0" << std::endl;
        assert(max > 0.0);
        return imax;
      }

#define SWAP(x,y) const int tmp=x; x=y; y=tmp;
      //#define SWAP(x,y)
      void lu_decompose(DenseMatrix & A, Vector<ComponentType> & b, std::vector<int> & r) const
      {
        // initialize row index array
        r.resize(rows);
        for(int i=0; i<rows; ++i)
          r[i]=i;

        //		std::cout << "r = ";
        //		for (int i=0; i<rows; i++)
        //			std::cout << r[i] << "  ";
        //		std::cout << std::endl;

        // row loop
        for(int k=0; k<rows-1; ++k)
          {
            // find pivot and swap indices
            int ipivot = pivot_index(A,k,k,r);
            SWAP(r[ipivot],r[k])

              //			std::cout << "r( "<< k <<") = ";
              //			for (int i=0; i<rows; i++)
              //				std::cout << r[i] << "  ";
              //			std::cout << std::endl;


              //			std::cout << std::endl;
              //			std::cout << "(A | b) (nach der " << k+1 << ".ten Zeilenvertauschung (oder auch nicht)) :" << std::endl;
              //			print( A, b, r );

              // elemination row loop
              for(int i=k+1; i<rows; ++i)
                {
                  ComponentType &a_ik = A(r[i],k);  // Wichtig: Diese Referenz "&" sorgt dafür, dass A(r[i],k) nach der Änderung von a_ik mit geändert wird!
                  //std::cout << "A(r[i],k) = " << A(r[i],k) << std::endl;
                  a_ik /= A(r[k],k);         // a_ik steht für "l_ik"(Literatur: Einträge von L) und wird im unteren Dreieck von A selbst gespeichert! Effizienz!
                  //std::cout << "A(r[i],k) = " << A(r[i],k) << std::endl;

                  // Jetzt kommt die Zeile-Subtraktion r[i]-te Zeile minus a_ik * r[k]-te Zeile:
                  // j ist hierbei der Spaltenlaufindex in der r[i]-ten Zeile
                  ComponentType *a_ij = &A(r[i],k); // Wichtig: Auch diese Referenzierung auf die Adresse von &A(r[i],k) sorgt dafür, dass A(r[i],k) nach der Änderung von a_ij mit geändert wird!
                  ComponentType *a_kj = &A(r[k],k);
                  for(int j=k+1; j<cols; ++j)
                    {
                      *(++a_ij) -=  a_ik * *(++a_kj);
                    }
                  b[r[i]] -= a_ik * b[r[k]];

                  //				std::cout << std::endl;
                  //				std::cout << "(A | b) (nach dem " << k+1 << "-ten Eliminationsschritt) :" << std::endl;
                  //				print( A, b, r );
                }
          }
      }

    public:
      void print(const DenseMatrix &A, const Vector<ComponentType> &b, const std::vector<int> &r) const
      {
        for(int ii=0; ii<rows; ii++)
          {
            std::cout << "( ";
            for(int kk=0; kk<cols; kk++)
              std::cout << std::setw(8) << A(r[ii], kk) << "  ";
            std::cout << " | ";
            std::cout << b[r[ii]] << "  ";
            std::cout << " )" << std::endl;
          }
      }

      // Constructor:
      DenseMatrix(const int _rows, const int _cols, const ComponentType def_val)
	: data( _rows*_cols, def_val )
	, rows(_rows)
        , cols(_cols)
        , row_sums(_rows, 0.0)
        , row_equilibrated( false )
      {}


      // Copy constructor:
      DenseMatrix( const DenseMatrix& other )
        : data( other.data )
        , rows( other.rows )
        , cols( other.cols )
        , row_sums( other.row_sums )
        , row_equilibrated( other.row_equilibrated )
      {}


      const int  & n_rows() const { return /**this.*/rows; }
      const int  & n_cols() const { return /**this.*/cols; }

      //int   n_rows(){ return rows; }
      //int   n_cols(){ return cols; }

      inline ComponentType & operator()(const int row, const int col)
      {
        assert(row < rows|| col < cols ||
               row >= 0 || col >= 0);

        return at(row,col);
      }



      inline const ComponentType & operator()(const int row, const int col) const
      {
        assert(row < rows|| col < cols ||
               row >= 0 || col >= 0);

        return at(row,col);
      }



      DenseMatrix& operator=( const DenseMatrix& other ){
        if( this != &other ){
          data.clear();
          data = other.data;
          rows = other.rows;
          cols = other.cols;
        }
        return *this;
      }


      Vector<ComponentType> operator* (const Vector<ComponentType> & x)
      {
        assert( int(x.size()) == this->rows );

        Vector<ComponentType> y( this->rows, 0 );
        for(int r=0; r<this->rows; ++r){
          for(int c=0; c<this->cols; ++c){
            y[r]+= at(r,c) * x[c];
          }
        }
        return y;
      }

      DenseMatrix operator* (const DenseMatrix & x)
      {
        assert(cols == x.n_rows());

        const int out_rows = rows;
        const int out_cols = x.n_cols();
        DenseMatrix y(out_rows, out_cols,0.0);
        for(int r=0; r<out_rows; ++r)
          for(int c=0; c<out_cols; ++c)
            for(int i=0; i<cols; ++i)
              y(r,c) += at(r,i) * x(i,c);

        return y;
      }

      void inverse(DenseMatrix<ComponentType> & X)
      {
        assert(cols == rows);
        assert(cols == X.n_cols());
        assert(rows == X.n_rows());

        for(int ii= 0 ;ii<cols;ii++){

          Vector<ComponentType> B(rows,0.0),X_row(rows,0.0);

          B[ii]=1.0; // canonical basis

          this->gauss_solver(X_row,B);
          for(int jj=0; jj<rows; jj++)
            X(jj,ii)=X_row[jj];
        }

      }

      void write_to_HDF5(const std::string & filename, const std::string & dataname){

        std::vector<UINT> dimensions(2,0);
        dimensions[0]=(UINT)rows;
        dimensions[1]=(UINT)cols;
        HDF5Tools::h5_Write( data,
                             filename,
                             dataname,
                             dimensions
                            );
      }

      void read_from_HDF5( const std::string & filename,
                           const std::string & dataname ){

        data.resize(0);
        rows=0;
        cols=0;

        Vector<UINT> local_offset(2,0);
        Vector<UINT> local_count(2,0);

        HDF5Tools::h5g_Read( data
                             , filename
                             , dataname
                             , local_offset
                             , local_count
                             );
        rows=local_count[0];
        cols=local_count[1];
      }


      void row_equilibration(){
        // This function must not be called twice!
        for(int i=0;i<this->rows;++i){
          row_sums[i] = 0.0;
          for(int j=0;j<this->cols;++j){
            row_sums[i] += std::abs( this->at(i,j) );
          }
          if( row_sums[i] < 1E-12 )
            std::cout << "BIG FAT WARNING: row sum almost zero" << std::endl;
          for(int j=0;j<this->cols;++j){
            this->at(i,j) /= row_sums[i];
          }
        }
        row_equilibrated = true;
      }


      void gauss_solver(Vector<ComponentType> & x, Vector<ComponentType> & b) const
      {

        if( row_equilibrated ){
          for(int i=0;i<this->rows;++i){
            b[i] /= row_sums[i];
          }
        }

        // make a copy of this matrix to hold the LU decomposition
        DenseMatrix LU(*this);
        // carve the solution out of the rhs
        x = b;

        // the row index array
        std::vector<int> r;

        // make LU decomposition
        lu_decompose(LU,x,r);

        // back insertion
        //		std::cout << std::endl;
        //		std::cout << "backward substitution:" << std::endl;
        for(int row=rows-1; row >= 0; --row)
          {
            //			std::cout << std::endl;
            //			std::cout << "Zeile " << row << ":" << std::endl;

            ComponentType *lu_rc = &LU(r[row],cols);

            //			std::cout << x[r[row]];

            for(int col=cols-1; col > row; --col)
              {
                //				Achtung!!! *(--lu_rc) zweimal zu benutzen fuer Debugging-Zwecke ist sehr gefaehrlich hier!!!
                //				std::cout << " - " << x[r[col]] << " * " << *(--lu_rc);
                x[r[row]] -= x[r[col]] * *(--lu_rc);

                // Korrektur 2:
                //				std::cout << " - " << x[r[col]] << " * " << LU(r[row],col);
                //				x[r[row]] -= x[r[col]] * LU(r[row],col);
              }
            //			std::cout << std::endl << " / " << LU(r[row],row) << "=" << std::endl;
            x[r[row]] /= LU(r[row],row);
            //			std::cout << x[r[row]] << " = " << "x[" << r[row] << "]" << std::endl;
          }

        // Korrektur 3:
        // Rücktransformation des Ergebnisvektors x:
        Vector<ComponentType> y(rows);
        for(int i=0; i<rows; i++)
          {
            y[i] = x[r[i]];
          }
        x = y;
      }
    };


    template<typename ComponentType>
    std::ostream & operator << ( std::ostream & os, const DenseMatrix<ComponentType> & x )
    {
      for(int r=0; r<x.n_rows(); ++r){
        os << "( ";
	for(int c=0; c<x.n_cols(); ++c)
          os << x(r,c) << " ";
        os << ") " << std::endl;
      }
      return os;
    }


  } // Gesis
} // Dune

#endif	/* DUNE_GESIS_DENSEMATRIX_HH */
