#include <complex>
#include <cmath>

class QuanObj
{
  private:
  
    size_t dim_i;  //  dimension of the operator;
    std::complex<double> *opr_i;  //  the quantum operator;
  
  public:
    
    QuanObj()
    {
      dim_i = 0;
    }
    
    QuanObj( size_t dim, std::complex<double> *opr )
    {
      dim_i = dim;
      opr_i = new std::complex<double> [ dim*dim ];
      for ( size_t i = 0; i < dim*dim; i++ )
        opr_i[i] = opr[i];
    }
    
    size_t get_dim()
    {
      return dim_i;
    }
    
    size_t set_opr( size_t dim, std::complex<double> *opr )
    {
      dim_i = dim;
      opr_i = new std::complex<double> [ dim*dim ];
      for ( size_t i = 0; i < dim*dim; i++ )
        opr_i[i] = opr[i];
      return 0;
    }
    
    std::complex<double> *get_opr()
    {
      std::complex<double> *opr;
      size_t dim = get_dim();
      opr = new std::complex<double> [ dim*dim ];
      for ( size_t i = 0; i < dim*dim; i++ )
        opr[i] = opr_i[i];
      return opr;
    }
};

class HamObj
{
  private:
  
    double  coeff;
    size_t  nbody;
    size_t  *pos_i;
    QuanObj *objlst;
    
  public:
    
    HamObj()
    {
      coeff = 0.0;
      nbody = 0;
    }
    
    HamObj( double coe )
    {
      coeff = coe;
      nbody = 0;
    }
    
    HamObj( size_t n, size_t *pos, QuanObj *objs )
    {
      coeff = 0.0;
      set_obj( n, pos, objs );
    }
    
    HamObj( double coe, size_t n, size_t *pos, QuanObj *objs )
    {
      coeff = coe;
      set_obj( n, pos, objs );
    }
    
    size_t set_coe( double coe )
    {
      coeff = coe;
      return 0;
    }
    
    double get_coe()
    {
      return coeff;
    }
    
    size_t get_nbody()
    {
      return nbody;
    }
    
    size_t set_obj( size_t n, size_t *pos, QuanObj *objs )
    {
      size_t i;
      nbody = n;
      pos_i = new size_t [n];
      objlst = new QuanObj [n];
      for ( i = 0; i < n; i++ )
      {
        pos_i[i] = pos[i];
        objlst[i] = objs[i];
      }
      return 0;
    }
    
    size_t *get_pos( size_t n, size_t *idx )
    {
      size_t i;
      size_t *pos;
      pos = new size_t [n];
      for ( i = 0; i < n; i++ )
        pos[i] = pos_i[ idx[i] ];
      return pos;
    }
    
    size_t *get_pos( size_t n )
    {
      size_t i, nex;
      size_t *pos;
      nex = std::min( n, nbody);
      pos = new size_t [ nex ];
      for ( i = 0; i < nex; i++ )
        pos[i] = pos_i[ i ];
      return pos;
    }
    
    QuanObj *get_obj( size_t n, size_t *idx )
    {
      size_t i;
      QuanObj *objs;
      
      objs = new QuanObj [n];
      for ( i = 0; i < n; i++ )
      {
        objs[i] = objlst[ idx[i] ];
      }
      
      return objs;
    }
    
    QuanObj *get_obj( size_t n )
    {
      size_t i, nex;
      QuanObj *objs;
      
      nex = std::min( n, nbody );
      objs = new QuanObj [ nex ];
      for ( i = 0; i < nex; i++ )
      {
        objs[i] = objlst[ i ];
      }
      
      return objs;
    }
    
    size_t append_obj( size_t n, size_t *pos, QuanObj *objs )
    {
      size_t i;
      size_t nnbody;
      size_t *pos_ex;
      QuanObj *obj_ex;
      
      nnbody = nbody + n;
      pos_ex = new size_t [nnbody];
      obj_ex = new QuanObj [nnbody];
      for ( i = 0; i < nbody; i++ )
      {
        pos_ex[ i ] = pos_i[ i ];
        obj_ex[ i ] = objlst[ i ];
      }
      for ( i = 0; i < n; i++ )
      {
        pos_ex[ nbody + i ] = pos[ i ];
        obj_ex[ nbody + i ] = objs[ i ];
      }
      
      nbody = nnbody;
      pos_i = pos_ex;
      objlst = obj_ex;
      
      return 0;
    }
};

class Ham
{
  private:
  
    size_t nspin;
    size_t nTerm;
    HamObj *hamlst;
    
  public:
    
    Ham()
    {
      nspin = 0;
      nTerm = 0;
    }
    
    Ham( size_t ns )
    {
      nspin = ns;
      nTerm = 0;
    }
    
    Ham( size_t ns, size_t nT, HamObj *hamobjs)
    {
      size_t i;
      nspin = ns;
      nTerm = nT;
      hamlst = new HamObj [nT];
      for ( i = 0; i < nT; i++ )
        hamlst[i] = hamobjs[i];
    }
    
    size_t set_nspin( size_t n )
    {
      nspin = n;
      return 0;
    }
    
    size_t get_nspin()
    {
      return nspin;
    }
    
    size_t get_nTerm()
    {
      return nTerm;
    }
    
    size_t set_obj( size_t nT, HamObj *obj )
    {
      size_t i;
      
      nTerm = nT;
      hamlst = new HamObj [nT];
      for ( i = 0; i < nT; i++ )
        hamlst[i] = obj[i];
      
      return 0;
    }
    
    size_t append_obj( size_t nT, HamObj *obj )
    {
      size_t i;
      size_t nnTerm;
      HamObj *hamlst_ex;
      
      nnTerm = nTerm + nT;
      hamlst_ex = new HamObj [nnTerm];
      for ( i = 0; i < nTerm; i++ )
        hamlst_ex[ i ] = hamlst[ i ];
      for ( i = 0; i < nT; i++ )
        hamlst_ex[ nTerm + i ] = obj[ i ];
      nTerm = nnTerm;
      hamlst = hamlst_ex;
      
      return 0;
    }
    
    HamObj *get_obj( size_t nT, size_t *idx )
    {
      size_t i;
      HamObj *hamobj;
      
      hamobj = new HamObj [nT];
      for ( i = 0; i < nT; i++ )
        hamobj[i] = hamlst[ idx[i] ];
      
      return hamobj;
    }
    
    size_t get_dim()
    {
      size_t     *nspin_lst, nspin_idx;
      size_t     nbody;
      size_t     *pos_i;
//      size_t     dim_i;
      QuanObj quanobj;
      HamObj  hamobj;
      
      size_t  dim;
      
      size_t     i, j, k, involved;
      
      nspin_lst = new size_t [nspin];
      for ( i = 0; i < nspin; i++ )
        nspin_lst[i] = 0;
      i = 0;
      hamobj = *get_obj( 1, &i );
      pos_i = hamobj.get_pos( 1 );
      nspin_lst[ 0 ] = pos_i[ 0 ];
      nspin_idx = 1;
      quanobj = *hamobj.get_obj( 1 );
      dim = quanobj.get_dim();
      
      for ( i = 0; i < nTerm; i++ )
      {
        hamobj = *get_obj( 1, &i );
        nbody  = hamobj.get_nbody();
        pos_i = hamobj.get_pos( nbody );
        for ( j = 0; j < nbody; j++ )
        {
          involved = 0;
          for ( k = 0; k < nspin_idx; k++ )
          {
            if ( pos_i[j] == nspin_lst[k] )
            {
              involved = 1;
              break;
            }
          }
          if ( involved == 0 )
          {
            if ( nspin_idx >= nspin )
            {
//              std::cout << "spin numbering exceed nspin!" << std::endl;
              return 0;
            }
            nspin_lst[ nspin_idx ] = pos_i[ j ];
            nspin_idx++;
            quanobj = *hamobj.get_obj( 1, &j );
            dim *= quanobj.get_dim();
          }
        }
      }
      
      return dim;
      
    }
};

