#ifndef MCNP2CAD_VOLUMES_H
#define MCNP2CAD_VOLUMES_H

#include <cstdlib>
#include "iGeom.h"
#include <armadillo>

// TODO: clean this igeom check function up
#define CHECK_BUF_SIZE 512
static char m_buf[CHECK_BUF_SIZE];

#define CHECK_IGEOM(err, msg) \
  do{/*std::cout << msg << std::endl;*/ if((err) != iBase_SUCCESS){     \
    std::cerr << "iGeom error (" << err << "): " << msg << std::endl;   \
    iGeom_getDescription( igm, m_buf, CHECK_BUF_SIZE);                  \
    std::cerr << " * " << m_buf << std::endl;                           \
     } } while(0) 

#endif /* MCNP2CAD_VOLUMES_H */


class Transform;
class SurfaceCard;

class SurfaceVolume{

protected:
  const Transform* transform;

public:
  SurfaceVolume( const Transform* transform_p = NULL):
    transform(transform_p)
  {}
  virtual ~SurfaceVolume(){}

  void setTransform( const Transform* transform_p ){ transform = transform_p; }
  
  virtual double getFarthestExtentFromOrigin( ) const = 0;
  virtual iBase_EntityHandle define( bool positive, iGeom_Instance& igm, double world_size );

protected:
  virtual iBase_EntityHandle getHandle( bool positive, iGeom_Instance& igm, double world_size ) = 0;
};


class VolumeCache;

/** 
 * Function to create an SurfaceVolume object from a SurfaceCard.
 * Created volumes are kept in a cache.  If the v parameter is null,
 * a default cache (static within volumes.cpp) will be used.  If multiple
 * MCNPInput instances are used in the course of a single program, then
 * there will need to be a way to get multiple VolumeCache objects;
 * that's not an issue right now, so I haven't coded it.
 */
extern 
SurfaceVolume& makeSurface( const SurfaceCard* card, VolumeCache* v = NULL );

extern 
iBase_EntityHandle makeWorldSphere( iGeom_Instance& igm, double world_size ); 

extern
iBase_EntityHandle applyTransform( const Transform& t, iGeom_Instance& igm, iBase_EntityHandle& e );

extern
iBase_EntityHandle applyReverseTransform( const Transform& tx, iGeom_Instance& igm, iBase_EntityHandle& e ) ;






class GeneralQuadraticSurface : public SurfaceVolume {

public:
  double A,B,C,D,E,F,G,H,J,K;

  GeneralQuadraticSurface( double _A, 
			   double _B,
			   double _C,
			   double _D,
			   double _E,
			   double _F,
			   double _G,
			   double _H,
			   double _J,
			   double _K ) :
    A(_A),B(_B),C(_C),D(_D),E(_E),F(_F),G(_G),H(_H),J(_J),K(_K)
  {     set_translation();
        set_rotation();
  }

protected:

  Vector3d translation;
  double rotation_mat[9];
  double extents[3];

  void set_translation();
  
  void set_rotation();

  virtual double getFarthestExtentFromOrigin ( ) const {
    return sqrt( extents[0]*extents[0] + extents[1]*extents[1] + extents[2]*extents[2]);
  }

virtual iBase_EntityHandle getHandle( bool positive, iGeom_Instance& igm, double world_size )
  { 

    std::cout << "Final Coeffs of GQ" << std::endl;

    std::cout << "A= " << A << std::endl;
    std::cout << "B= " << B << std::endl;
    std::cout << "C= " << C << std::endl;
    std::cout << "D= " << D << std::endl;
    std::cout << "E= " << E << std::endl;
    std::cout << "F= " << F << std::endl;
    std::cout << "K= " << K << std::endl;


    iBase_EntityHandle gq_handle;
    int igm_result=0;
    iGeom_GQ(igm,A,B,C,D,E,F,G,H,J,K,world_size,&gq_handle,&igm_result);
    CHECK_IGEOM( igm_result, "Creating intial GQ");

    Transform gq_transform(rotation_mat,translation);

    //move back to original orientation
    applyReverseTransform( gq_transform, igm, gq_handle);

    double xmin,xmax,ymin,ymax,zmin,zmax;
    iGeom_getBoundBox(igm,&xmin,&ymin,&zmin,&xmax,&ymax,&zmax,&igm_result);
    extents[0] = ( fabs(xmin) > fabs(xmax) ) ? xmin : xmax;
    extents[1] = ( fabs(ymin) > fabs(ymax) ) ? ymin : ymax;
    extents[2] = ( fabs(zmin) > fabs(zmax) ) ? zmin : zmax;

    return gq_handle;
  }




};
