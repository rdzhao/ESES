#include "MarchingCubes.h"

int MarchingCubes::read_bounding_box(std::string filename)
{
  std::fstream file;
  file.open(filename.c_str(), std::ios::in);
  double xl, yf, zb;
  double xr, yb, zt;

  file>>xl>>yf>>zb;
  file>>xr>>yb>>zt;
  file>>_x_size>>_y_size>>_z_size;
  
  std::vector<double> tmp_v;
  tmp_v.push_back(xl);
  tmp_v.push_back(yf);
  tmp_v.push_back(zb);
  _boun_box.push_back(tmp_v);
  tmp_v.clear();
  tmp_v.push_back(xr);
  tmp_v.push_back(yb);
  tmp_v.push_back(zt);
  _boun_box.push_back(tmp_v);
  
  return 1;
}

int MarchingCubes::read_grid_info(std::string filename)
{
  std::fstream file;
  file.open(filename.c_str(), std::ios::in);
  
  std::vector<double> tmp_v;
  int tmp_val;
  for(int k=0; k<_z_size; k++)
    for(int j=0; j<_y_size; j++)
      for(int i=0; i<_x_size; i++)
	{
	  file>>tmp_val;
	  tmp_v.push_back(tmp_val);
	  file>>tmp_val;
	  tmp_v.push_back(tmp_val);
	  file>>tmp_val;
	  tmp_v.push_back(tmp_val);
	  file>>tmp_val;
	  tmp_v.push_back(tmp_val);
	  _grid_info.push_back(tmp_v);
	  tmp_v.clear();
	}

  return 1;
}

int MarchingCubes::read_intersection_info(std::string filename)
{
  std::fstream file;
  file.open(filename.c_str(), std::ios::in);
  
  std::vector<int> tmp_i;
  std::vector<double> tmp_d;
  std::string line;
  
  while(getline(file, line)){
      std::stringstream ss(line);
      int token_i;
      double token_d;
      
      for(int i=0; i<6; i++){
	  ss>>token_i;
	  tmp_i.push_back(token_i);
	}
      _inte_dual.push_back(tmp_i);
      tmp_i.clear();
      
      for(int i=0; i<3; i++){
	  ss>>token_d;
	  tmp_d.push_back(token_d);
	}
      _inte_verts.push_back(tmp_d);
      tmp_d.clear();
      
      for(int i=0; i<3; i++){
	  ss>>token_d;
	  tmp_d.push_back(token_d);
	}
      _inte_normal.push_back(tmp_d);
      tmp_d.clear();
    }

  return 1;
}

int MarchingCubes::process(){
    set_bounding_box();
    set_grid();
    set_intersection();

    initialize_verts();
    save_intersection_info();
    run();
    write_OBJ_file();

    return 1;
}

int MarchingCubes::set_bounding_box(){
    std::vector<double> tmp;
    tmp = {_pw->get_basic_info().t_min_x, _pw->get_basic_info().t_min_y, _pw->get_basic_info().t_min_z};
    _boun_box.push_back(tmp);
    tmp = {_pw->get_basic_info().t_max_x, _pw->get_basic_info().t_max_y, _pw->get_basic_info().t_max_z};
    _boun_box.push_back(tmp);

    _x_size = _pw->get_basic_info().t_x_num;
    _y_size = _pw->get_basic_info().t_y_num;
    _z_size = _pw->get_basic_info().t_z_num;

    return 1;
}

int MarchingCubes::set_grid(){
    for(int i=0; i<_pw->get_grid_status().size(); ++i){
        _data.push_back(_pw->get_grid_status()[i] ? 1 : -1);
    }
    return 1;
}

int MarchingCubes::set_intersection(){
    std::vector<int> tmp_int;
    std::vector<double> tmp_double;

    for(int i=0; i<_pw->get_intersection_info().m_intersection.size(); ++i){
        CVector3d& gv_in = _pw->get_intersection_info().m_grid_inside[i];
        CVector3d& gv_out = _pw->get_intersection_info().m_grid_outside[i];
        tmp_int = { static_cast<int>(gv_in.x()),
                    static_cast<int>(gv_in.y()),
                    static_cast<int>(gv_in.z()),
                    static_cast<int>(gv_out.x()),
                    static_cast<int>(gv_out.y()),
                    static_cast<int>(gv_out.z())};
        _inte_dual.push_back(tmp_int);

        CVector3d& pos = _pw->get_intersection_info().m_intersection[i];
        CVector3d& normal = _pw->get_intersection_info().m_normal[i];
        _inte_verts.push_back({pos.x(), pos.y(), pos.z()});
        _inte_normal.push_back({normal.x(), normal.y(), normal.z()});
    }
    return 1;
}

int MarchingCubes::save_grid_info()
{
  for(int i=0; i<_grid_info.size(); i++)
    _data.push_back(_grid_info[i][3]);

  return 1;
}

int MarchingCubes::initialize_verts()
{
  for(int i=0; i<_x_size*_y_size*_z_size; i++)
    _x_verts.push_back(-1);
  for(int i=0; i<_x_size*_y_size*_z_size; i++)
    _y_verts.push_back(-1);
  for(int i=0; i<_x_size*_y_size*_z_size; i++)
    _z_verts.push_back(-1);

  return 1;
}

int MarchingCubes::add_x_vert(const int i)
{
  _nvertices++;
  //std::cout<<"before entering add_x"<<std::endl;
  //std::cout<<"x: "<<i<<std::endl;
  Vertex v;
  v.x=_inte_verts[i][0];
  v.y=_inte_verts[i][1];
  v.z=_inte_verts[i][2];
  v.nx=_inte_normal[i][0];
  v.ny=_inte_normal[i][1];
  v.nz=_inte_normal[i][2];
  //std::cout<<_vertices.size()<<std::endl;
  
  _vertices.push_back(v);
  //std::cout<<"x: "<<i<<std::endl;
  return _nvertices-1;
}

int MarchingCubes::add_y_vert(const int i)
{
  _nvertices++;
  
  Vertex v;
  v.x=_inte_verts[i][0];
  v.y=_inte_verts[i][1];
  v.z=_inte_verts[i][2];
  v.nx=_inte_normal[i][0];
  v.ny=_inte_normal[i][1];
  v.nz=_inte_normal[i][2];
  _vertices.push_back(v);
  //std::cout<<"y: "<<i<<std::endl;
  return _nvertices-1;
}

int MarchingCubes::add_z_vert(const int i)
{
  _nvertices++;
  
  Vertex v;
  v.x=_inte_verts[i][0];
  v.y=_inte_verts[i][1];
  v.z=_inte_verts[i][2];
  v.nx=_inte_normal[i][0];
  v.ny=_inte_normal[i][1];
  v.nz=_inte_normal[i][2];
  _vertices.push_back(v);
  //std::cout<<"z: "<<i<<std::endl;
  return _nvertices-1;
}

int MarchingCubes::add_c_vert()
{
  //std::cout<<"in add_c"<<std::endl;
  _nvertices++;
  Vertex v;
  v.x=0;
  v.y=0;
  v.z=0;
  v.nx=0;
  v.ny=0;
  v.nz=0;
  int u=0;
  int vid;

  vid=get_x_vert(_i, _j, _k);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_y_vert(_i+1, _j, _k);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_x_vert(_i, _j+1, _k);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_y_vert(_i, _j, _k);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_x_vert(_i, _j, _k+1);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_y_vert(_i+1, _j, _k+1);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_x_vert(_i, _j+1, _k+1);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_y_vert(_i, _j, _k+1);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_z_vert(_i, _j, _k);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_z_vert(_i+1, _j, _k);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_z_vert(_i+1, _j+1, _k);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }
  vid=get_z_vert(_i, _j+1, _k);
  if(vid!=-1)
    {
      u++;
      const Vertex vert=_vertices[vid];
      v.x+=vert.x;
      v.y+=vert.y;
      v.z+=vert.z;
      v.nx+=vert.nx;
      v.ny+=vert.ny;
      v.nz+=vert.nz;
    }

  v.x/=u;
  v.y/=u;
  v.z/=u;
  //std::cout<<u<<std::endl;
  double sr=sqrt(v.nx*v.nx+v.ny*v.ny+v.nz*v.nz);
  v.nx/=sr;
  v.ny/=sr;
  v.nz/=sr;
  _vertices.push_back(v);
  
  return _nvertices-1;
}

int MarchingCubes::set_x_vert(const int val, const int i, const int j, const int k)
{
  //std::cout<<"before"<<std::endl;
  //std::cout<<i+j*_x_size+k*_x_size*_y_size<<std::endl;
  _x_verts[i+j*_x_size+k*_x_size*_y_size]=val;
  //std::cout<<val<<std::endl;
  //std::cout<<"after"<<std::endl;
  return 1;
}

int MarchingCubes::set_y_vert(const int val, const int i, const int j, const int k)
{
  _y_verts[i+j*_x_size+k*_x_size*_y_size]=val;
  return 1;
}

int MarchingCubes::set_z_vert(const int val, const int i, const int j, const int k)
{
  _z_verts[i+j*_x_size+k*_x_size*_y_size]=val;
  return 1;
}

int MarchingCubes::save_intersection_info()
{
  for(int i=0; i<_inte_dual.size(); i++)
    {
      //std::cout<<"into for"<<i<<std::endl;
      if(_inte_dual[i][0]!=_inte_dual[i][3])
	{
	  //std::cout<<"before saving"<<std::endl;
	  if(_inte_dual[i][0]>_inte_dual[i][3])
	    {	
	      int tmp=add_x_vert(i);
	      set_x_vert(tmp, _inte_dual[i][3], _inte_dual[i][4], _inte_dual[i][5]);
	    }
	  else
	    {
	      int tmp=add_x_vert(i);
	      set_x_vert(tmp, _inte_dual[i][0], _inte_dual[i][1], _inte_dual[i][2]);
	    }
	  //std::cout<<"after saving"<<std::endl;
	}
      //std::cout<<"out for"<<std::endl;
      
      else if(_inte_dual[i][1]!=_inte_dual[i][4])
	{
	  if(_inte_dual[i][1]>_inte_dual[i][4])	
	    set_y_vert(add_y_vert(i), _inte_dual[i][3], _inte_dual[i][4], _inte_dual[i][5]);
	  else
	    set_y_vert(add_y_vert(i), _inte_dual[i][0], _inte_dual[i][1], _inte_dual[i][2]);
	}
      else if(_inte_dual[i][2]!=_inte_dual[i][5])
	{
	  if(_inte_dual[i][2]>_inte_dual[i][5])	
	    set_z_vert(add_z_vert(i), _inte_dual[i][3], _inte_dual[i][4], _inte_dual[i][5]);
	  else
	    set_z_vert(add_z_vert(i), _inte_dual[i][0], _inte_dual[i][1], _inte_dual[i][2]);
	}
      
      //std::cout<<"dual index: "<<i<<std::endl;
    }
  //std::cout<<"out save inte"<<std::endl;
  return 1;
}

int MarchingCubes::test_face(signed char face)
{
  double A,B,C,D ;

  switch( face )
    {
    case -1 : case 1 :  A = _cube[0] ;  B = _cube[4] ;  C = _cube[5] ;  D = _cube[1] ;  break ;
    case -2 : case 2 :  A = _cube[1] ;  B = _cube[5] ;  C = _cube[6] ;  D = _cube[2] ;  break ;
    case -3 : case 3 :  A = _cube[2] ;  B = _cube[6] ;  C = _cube[7] ;  D = _cube[3] ;  break ;
    case -4 : case 4 :  A = _cube[3] ;  B = _cube[7] ;  C = _cube[4] ;  D = _cube[0] ;  break ;
    case -5 : case 5 :  A = _cube[0] ;  B = _cube[3] ;  C = _cube[2] ;  D = _cube[1] ;  break ;
    case -6 : case 6 :  A = _cube[4] ;  B = _cube[7] ;  C = _cube[6] ;  D = _cube[5] ;  break ;
    default : printf( "Invalid face code %d\n", face ) ;  A = B = C = D = 0 ;
    };

  if( fabs( A*C - B*D ) < FLT_EPSILON )
    return face >= 0 ;
  return face * A * ( A*C - B*D ) >= 0  ;  // face and A invert signs
}

bool MarchingCubes::test_interior(signed char s)
{
  double t, At=0, Bt=0, Ct=0, Dt=0, a, b ;
  char  test =  0 ;
  char  edge = -1 ; // reference edge of the triangulation

  switch( _case )
    {
    case  4 :
    case 10 :
      a = ( _cube[4] - _cube[0] ) * ( _cube[6] - _cube[2] ) - ( _cube[7] - _cube[3] ) * ( _cube[5] - _cube[1] ) ;
      b =  _cube[2] * ( _cube[4] - _cube[0] ) + _cube[0] * ( _cube[6] - _cube[2] )
	- _cube[1] * ( _cube[7] - _cube[3] ) - _cube[3] * ( _cube[5] - _cube[1] ) ;
      t = - b / (2*a) ;
      if( t<0 || t>1 ) return s>0 ;

      At = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
      Bt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
      Ct = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
      Dt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
      break ;

    case  6 :
    case  7 :
    case 12 :
    case 13 :
      switch( _case )
	{
	case  6 : edge = test6 [_config][2] ; break ;
	case  7 : edge = test7 [_config][4] ; break ;
	case 12 : edge = test12[_config][3] ; break ;
	case 13 : edge = tiling13_5_1[_config][_subconfig][0] ; break ;
	}
      switch( edge )
	{
	case  0 :
	  t  = _cube[0] / ( _cube[0] - _cube[1] ) ;
	  At = 0 ;
	  Bt = _cube[3] + ( _cube[2] - _cube[3] ) * t ;
	  Ct = _cube[7] + ( _cube[6] - _cube[7] ) * t ;
	  Dt = _cube[4] + ( _cube[5] - _cube[4] ) * t ;
	  break ;
	case  1 :
	  t  = _cube[1] / ( _cube[1] - _cube[2] ) ;
	  At = 0 ;
	  Bt = _cube[0] + ( _cube[3] - _cube[0] ) * t ;
	  Ct = _cube[4] + ( _cube[7] - _cube[4] ) * t ;
	  Dt = _cube[5] + ( _cube[6] - _cube[5] ) * t ;
	  break ;
	case  2 :
	  t  = _cube[2] / ( _cube[2] - _cube[3] ) ;
	  At = 0 ;
	  Bt = _cube[1] + ( _cube[0] - _cube[1] ) * t ;
	  Ct = _cube[5] + ( _cube[4] - _cube[5] ) * t ;
	  Dt = _cube[6] + ( _cube[7] - _cube[6] ) * t ;
	  break ;
	case  3 :
	  t  = _cube[3] / ( _cube[3] - _cube[0] ) ;
	  At = 0 ;
	  Bt = _cube[2] + ( _cube[1] - _cube[2] ) * t ;
	  Ct = _cube[6] + ( _cube[5] - _cube[6] ) * t ;
	  Dt = _cube[7] + ( _cube[4] - _cube[7] ) * t ;
	  break ;
	case  4 :
	  t  = _cube[4] / ( _cube[4] - _cube[5] ) ;
	  At = 0 ;
	  Bt = _cube[7] + ( _cube[6] - _cube[7] ) * t ;
	  Ct = _cube[3] + ( _cube[2] - _cube[3] ) * t ;
	  Dt = _cube[0] + ( _cube[1] - _cube[0] ) * t ;
	  break ;
	case  5 :
	  t  = _cube[5] / ( _cube[5] - _cube[6] ) ;
	  At = 0 ;
	  Bt = _cube[4] + ( _cube[7] - _cube[4] ) * t ;
	  Ct = _cube[0] + ( _cube[3] - _cube[0] ) * t ;
	  Dt = _cube[1] + ( _cube[2] - _cube[1] ) * t ;
	  break ;
	case  6 :
	  t  = _cube[6] / ( _cube[6] - _cube[7] ) ;
	  At = 0 ;
	  Bt = _cube[5] + ( _cube[4] - _cube[5] ) * t ;
	  Ct = _cube[1] + ( _cube[0] - _cube[1] ) * t ;
	  Dt = _cube[2] + ( _cube[3] - _cube[2] ) * t ;
	  break ;
	case  7 :
	  t  = _cube[7] / ( _cube[7] - _cube[4] ) ;
	  At = 0 ;
	  Bt = _cube[6] + ( _cube[5] - _cube[6] ) * t ;
	  Ct = _cube[2] + ( _cube[1] - _cube[2] ) * t ;
	  Dt = _cube[3] + ( _cube[0] - _cube[3] ) * t ;
	  break ;
	case  8 :
	  t  = _cube[0] / ( _cube[0] - _cube[4] ) ;
	  At = 0 ;
	  Bt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
	  Ct = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
	  Dt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
	  break ;
	case  9 :
	  t  = _cube[1] / ( _cube[1] - _cube[5] ) ;
	  At = 0 ;
	  Bt = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
	  Ct = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
	  Dt = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
	  break ;
	case 10 :
	  t  = _cube[2] / ( _cube[2] - _cube[6] ) ;
	  At = 0 ;
	  Bt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
	  Ct = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
	  Dt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
	  break ;
	case 11 :
	  t  = _cube[3] / ( _cube[3] - _cube[7] ) ;
	  At = 0 ;
	  Bt = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
	  Ct = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
	  Dt = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
	  break ;
	default : printf( "Invalid edge %d\n", edge ) ;  break ;
	}
      break ;

    default : printf( "Invalid ambiguous case %d\n", _case ) ;  break ;
    }

  if( At >= 0 ) test ++ ;
  if( Bt >= 0 ) test += 2 ;
  if( Ct >= 0 ) test += 4 ;
  if( Dt >= 0 ) test += 8 ;
  switch( test )
    {
    case  0 : return s>0 ;
    case  1 : return s>0 ;
    case  2 : return s>0 ;
    case  3 : return s>0 ;
    case  4 : return s>0 ;
    case  5 : if( At * Ct - Bt * Dt <  FLT_EPSILON ) return s>0 ; break ;
    case  6 : return s>0 ;
    case  7 : return s<0 ;
    case  8 : return s>0 ;
    case  9 : return s>0 ;
    case 10 : if( At * Ct - Bt * Dt >= FLT_EPSILON ) return s>0 ; break ;
    case 11 : return s<0 ;
    case 12 : return s>0 ;
    case 13 : return s<0 ;
    case 14 : return s<0 ;
    case 15 : return s<0 ;
    }

  return s<0 ;
}

int MarchingCubes::add_triangle(const char* trig, int n, int v12)
{
  std::vector<int> tv;
  tv.resize(3);

  for(int i=0; i<3*n; i++)
    {
      switch(trig[i])
	{
	case 0: tv[i%3]=get_x_vert(_i, _j, _k); break;
	case 1: tv[i%3]=get_y_vert(_i+1, _j, _k); break;
	case 2: tv[i%3]=get_x_vert(_i, _j+1, _k); break;
	case 3: tv[i%3]=get_y_vert(_i, _j, _k); break;
	case 4: tv[i%3]=get_x_vert(_i, _j, _k+1); break;
	case 5: tv[i%3]=get_y_vert(_i+1, _j, _k+1); break;
	case 6: tv[i%3]=get_x_vert(_i, _j+1, _k+1); break;
	case 7: tv[i%3]=get_y_vert(_i, _j, _k+1); break;
	case 8: tv[i%3]=get_z_vert(_i, _j, _k); break;
	case 9: tv[i%3]=get_z_vert(_i+1, _j, _k); break;
	case 10: tv[i%3]=get_z_vert(_i+1, _j+1, _k); break;
	case 11: tv[i%3]=get_z_vert(_i, _j+1, _k); break;
	case 12: tv[i%3]=v12; break;
	default: break;
	};

      //if(tv[i%3]==-1)
      //std::cout<<"Invalid triangle "<<_ntriangles+1<<std::endl;
      
      if(i%3==2)
	{
	  _ntriangles++;
	  Triangle t;
	  t.v1=tv[0];
	  t.v2=tv[1];
	  t.v3=tv[2];
	  _triangles.push_back(t);
	}
    }

  return 1;
}

int MarchingCubes::process_cube()
{
  int v12=-1 ;
  _case=cases[_lut_entry][0];
  _config=cases[_lut_entry][1];
  _subconfig=0;

  switch(_case)
    {
    case 0:
      break;
      
    case 1:
      add_triangle(tiling1[_config], 1);
      break;

    case 2:
      add_triangle(tiling2[_config], 2);
      break;

    case 3:
      if(test_face(test3[_config]))
	add_triangle(tiling3_2[_config], 4); // 3.2
      else
	add_triangle(tiling3_1[_config], 2); // 3.1
      break;

    case 4:
      if(test_interior(test4[_config]))
	add_triangle(tiling4_1[_config], 2); // 4.1.1
      else
	add_triangle(tiling4_2[_config], 6); // 4.1.2
      break;

    case 5:
      add_triangle(tiling5[_config], 3);
      break;

    case 6:
      if(test_face(test6[_config][0]))
	add_triangle(tiling6_2[_config], 5); // 6.2
      else
	{
	  if(test_interior(test6[_config][1]))
	    add_triangle(tiling6_1_1[_config], 3); // 6.1.1
	  else
	    {
	      v12=add_c_vert();
	      add_triangle(tiling6_1_2[_config], 9, v12); // 6.1.2
	    }
	}
      break;

    case 7:
      if(test_face(test7[_config][0])) _subconfig+=1;
      if(test_face(test7[_config][1])) _subconfig+=2;
      if(test_face(test7[_config][2])) _subconfig+=4 ;
      switch(_subconfig)
	{
	case 0:
	  add_triangle(tiling7_1[_config], 3); break;
	case 1:
	  add_triangle(tiling7_2[_config][0], 5); break;
	case 2:
	  add_triangle(tiling7_2[_config][1], 5); break;
	case 3:
	  v12=add_c_vert();
	  add_triangle(tiling7_3[_config][0], 9, v12); break;
	case 4:
	  add_triangle(tiling7_2[_config][2], 5); break;
	case 5:
	  v12=add_c_vert();
	  add_triangle(tiling7_3[_config][1], 9, v12); break;
	case 6:
	  v12=add_c_vert();
	  add_triangle(tiling7_3[_config][2], 9, v12); break;
	case 7:
	  if(test_interior(test7[_config][3]))
	    add_triangle(tiling7_4_2[_config], 9);
	  else
	    add_triangle(tiling7_4_1[_config], 5);
	  break;
	};
      break;

    case 8:
      add_triangle(tiling8[_config], 2);
      break;

    case 9:
      add_triangle(tiling9[_config], 4);
      break;

    case 10:
      if(test_face(test10[_config][0]))
	{
	  if(test_face(test10[_config][1]))
	    add_triangle(tiling10_1_1_[_config], 4); // 10.1.1
	  else
	    {
	      v12=add_c_vert();
	      add_triangle(tiling10_2[_config], 8, v12); // 10.2
	    }
	}
      else
	{
	  if(test_face(test10[_config][1]))
	    {
	      v12=add_c_vert();
	      add_triangle(tiling10_2_[_config], 8, v12); // 10.2
	    }
	  else
	    {
	      if(test_interior(test10[_config][2]))
		add_triangle(tiling10_1_1[_config], 4); // 10.1.1
	      else
		add_triangle(tiling10_1_2[_config], 8); // 10.1.2
	    }
	}
      break;

    case 11:
      add_triangle(tiling11[_config], 4);
      break;

    case 12:
      if(test_face(test12[_config][0]))
	{
	  if(test_face(test12[_config][1]))
	    add_triangle(tiling12_1_1_[_config], 4); // 12.1.1
	  else
	    {
	      v12=add_c_vert();
	      add_triangle(tiling12_2[_config], 8, v12); // 12.2
	    }
	}
      else
	{
	  if(test_face(test12[_config][1]))
	    {
	      v12=add_c_vert();
	      add_triangle(tiling12_2_[_config], 8, v12); // 12.2
	    }
	  else
	    {
	      if(test_interior(test12[_config][2]))
		add_triangle(tiling12_1_1[_config], 4); // 12.1.1
	      else
		add_triangle(tiling12_1_2[_config], 8); // 12.1.2
	    }
	}
      break;

    case 13:
      if(test_face(test13[_config][0])) _subconfig+=1;
      if(test_face(test13[_config][1])) _subconfig+=2;
      if(test_face(test13[_config][2])) _subconfig+=4;
      if(test_face(test13[_config][3])) _subconfig+=8;
      if(test_face(test13[_config][4])) _subconfig+=16;
      if(test_face(test13[_config][5])) _subconfig+=32;
      switch(subconfig13[_subconfig])
	{
	case 0:/* 13.1 */
	  add_triangle(tiling13_1[_config], 4); break;

	case 1:/* 13.2 */
	  add_triangle(tiling13_2[_config][0], 6); break;
	case 2:/* 13.2 */
	  add_triangle(tiling13_2[_config][1], 6); break;
	case 3:/* 13.2 */
	  add_triangle(tiling13_2[_config][2], 6); break;
	case 4:/* 13.2 */
	  add_triangle(tiling13_2[_config][3], 6); break;
	case 5:/* 13.2 */
	  add_triangle(tiling13_2[_config][4], 6); break;
	case 6:/* 13.2 */
	  add_triangle(tiling13_2[_config][5], 6); break;

	case 7:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][0], 10, v12); break;
	case 8:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][1], 10, v12); break;
	case 9:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][2], 10, v12); break;
	case 10:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][3], 10, v12); break;
	case 11:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][4], 10, v12); break;
	case 12:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][5], 10, v12); break;
	case 13:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][6], 10, v12); break;
	case 14:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][7], 10, v12); break;
	case 15:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][8], 10, v12); break;
	case 16:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][9], 10, v12); break;
	case 17:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][10], 10, v12); break;
	case 18:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3[_config][11], 10, v12); break;

	case 19:/* 13.4 */
	  v12=add_c_vert();
	  add_triangle(tiling13_4[_config][0], 12, v12); break;
	case 20:/* 13.4 */
	  v12=add_c_vert();
	  add_triangle(tiling13_4[_config][1], 12, v12); break;
	case 21:/* 13.4 */
	  v12=add_c_vert();
	  add_triangle(tiling13_4[_config][2], 12, v12); break;
	case 22:/* 13.4 */
	  v12=add_c_vert();
	  add_triangle(tiling13_4[_config][3], 12, v12); break;

	case 23:/* 13.5 */
	  _subconfig=0;
	  if(test_interior(test13[_config][6]))
	    add_triangle(tiling13_5_1[_config][0], 6);
	  else
	    add_triangle(tiling13_5_2[_config][0], 10);
	  break;
	case 24:/* 13.5 */
	  _subconfig=1;
	  if(test_interior(test13[_config][6]))
	    add_triangle(tiling13_5_1[_config][1], 6);
	  else
	    add_triangle(tiling13_5_2[_config][1], 10);
	  break;
	case 25:/* 13.5 */
	  _subconfig=2;
	  if(test_interior(test13[_config][6]))
	    add_triangle(tiling13_5_1[_config][2], 6);
	  else
	    add_triangle(tiling13_5_2[_config][2], 10);
	  break;
	case 26:/* 13.5 */
	  _subconfig=3;
	  if(test_interior(test13[_config][6]))
	    add_triangle(tiling13_5_1[_config][3], 6);
	  else
	    add_triangle(tiling13_5_2[_config][3], 10);
	  break;

	case 27:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][0], 10, v12); break;
	case 28:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][1], 10, v12); break;
	case 29:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][2], 10, v12); break;
	case 30:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][3], 10, v12); break;
	case 31:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][4], 10, v12); break;
	case 32:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][5], 10, v12); break;
	case 33:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][6], 10, v12); break;
	case 34:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][7], 10, v12); break;
	case 35:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][8], 10, v12); break;
	case 36:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][9], 10, v12); break;
	case 37:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][10], 10, v12); break;
	case 38:/* 13.3 */
	  v12=add_c_vert();
	  add_triangle(tiling13_3_[_config][11], 10, v12); break;

	case 39:/* 13.2 */
	  add_triangle(tiling13_2_[_config][0], 6); break;
	case 40:/* 13.2 */
	  add_triangle(tiling13_2_[_config][1], 6); break;
	case 41:/* 13.2 */
	  add_triangle(tiling13_2_[_config][2], 6); break;
	case 42:/* 13.2 */
	  add_triangle(tiling13_2_[_config][3], 6); break;
	case 43:/* 13.2 */
	  add_triangle(tiling13_2_[_config][4], 6); break;
	case 44:/* 13.2 */
	  add_triangle(tiling13_2_[_config][5], 6); break;

	case 45:/* 13.1 */
	  add_triangle(tiling13_1_[_config], 4); break;

	default :
	  std::cout<<"Marching Cubes: Impossible case 13?"<<std::endl;
	}
      break;

    case 14:
      add_triangle(tiling14[_config], 4);
      break;
    }
  return 1;
}

int MarchingCubes::run()
{
  for( _k = 0 ; _k < _z_size-1 ; _k++ )
    for( _j = 0 ; _j < _y_size-1 ; _j++ )
      for( _i = 0 ; _i < _x_size-1 ; _i++ )
	{
	  _lut_entry = 0 ;
	  for( int p = 0 ; p < 8 ; ++p )
	    {
	      _cube[p] = get_data( _i+((p^(p>>1))&1), _j+((p>>1)&1), _k+((p>>2)&1) ) ;
	      if( fabs( _cube[p] ) < FLT_EPSILON ) _cube[p] = FLT_EPSILON ;
	      if( _cube[p] > 0 ) _lut_entry += 1 << p ;
	    }
	  process_cube( ) ;
	}
  return 1;
}

int MarchingCubes::triangle_normal()
{
	for (int i = 0; i < _triangles.size(); i++)
	{
	  //std::cout<<i<<std::endl;
		std::vector<double> v1;
		v1.push_back(_vertices[_triangles[i].v3].x - _vertices[_triangles[i].v1].x);
		v1.push_back(_vertices[_triangles[i].v3].y - _vertices[_triangles[i].v1].y);
		v1.push_back(_vertices[_triangles[i].v3].z - _vertices[_triangles[i].v1].z);
		//std::cout<<i<<std::endl;
		std::vector<double> v2;
		v2.push_back(_vertices[_triangles[i].v2].x - _vertices[_triangles[i].v3].x);
		v2.push_back(_vertices[_triangles[i].v2].y - _vertices[_triangles[i].v3].y);
		v2.push_back(_vertices[_triangles[i].v2].z - _vertices[_triangles[i].v3].z);
		std::vector<double> cross;
		cross.push_back(v1[1] * v2[2] - v1[2] * v2[1]);
		cross.push_back(-v1[0] * v2[2] + v1[2] * v2[0]);
		cross.push_back(v1[0] * v2[1] - v1[1] * v2[0]);
		double norm = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
		cross[0] /= norm;
		cross[1] /= norm;
		cross[2] /= norm;
		_triangle_normal.push_back(cross);
		//std::cout<<i<<std::endl;
	}

	return 1;
}

int MarchingCubes::write_OBJ_file()
{
  std::fstream file;
  file.open("marching_cubes.obj", std::ios::out);
  std::cout<<"in func write, vertices size: "<<_vertices.size()<<std::endl;
  for(int i=0; i<_vertices.size(); i++)
    file<<"v "
	<<_vertices[i].x<<" "
	<<_vertices[i].y<<" "
	<<_vertices[i].z<<std::endl;
  //for(int i=0; i<_vertices.size(); i++)
  //  file<<"vn "
	//<<_vertices[i].nx<<" "
	//<<_vertices[i].ny<<" "
	//<<_vertices[i].nz<<std::endl;
  for(int i=0; i<_triangles.size(); i++)
    file<<"f "
	<<_triangles[i].v1+1<<" "
	<<_triangles[i].v3+1<<" "
	<<_triangles[i].v2+1<<std::endl;
  
  file.close();

  return 1;
}

int MarchingCubes::write_OFF_file()
{
	std::fstream file;
  	file.open("MC_result.off", std::ios::out);
	std::cout<<"in func write, vertices size: "<<_vertices.size()<<std::endl;
  	
	file<<"OFF"<<std::endl;
	file<<_vertices.size()<<" "<<_triangles.size()<<" 0"<<std::endl;  
	for(int i=0; i<_vertices.size(); i++)
    	file<<_vertices[i].x<<" "
			<<_vertices[i].y<<" "
			<<_vertices[i].z<<std::endl;

  	for(int i=0; i<_triangles.size(); i++)
    	file<<"3 "
			<<_triangles[i].v1<<" "
			<<_triangles[i].v3<<" "
			<<_triangles[i].v2<<std::endl;
  
  	file.close();

  	return 1;
}

int MarchingCubes::write_VTK_file()
{
	std::fstream file;
	file.open("MC_result.vtk", std::ios::out);
	std::cout<<"in func write VTK, vertices size: "<<_vertices.size()<<std::endl;

	file<<"# vtk DataFile Version 3.0"<<std::endl;
	file<<"MC_result"<<std::endl;
	file<<"ASCII"<<std::endl;
	file<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	file<<std::endl;

	file<<"POINTS "<<_vertices.size()<<" double"<<std::endl;
	for(int i=0; i<_vertices.size(); ++i){
		file<<_vertices[i].x<<" "
			<<_vertices[i].y<<" "
			<<_vertices[i].z<<std::endl;
	}
	file<<std::endl;

	file<<"CELLS "<<_triangles.size()<<" "<<4*_triangles.size()<<std::endl;
	for(int i=0; i<_triangles.size(); ++i){
		file<<"3 "
			<<_triangles[i].v1<<" "
			<<_triangles[i].v2<<" "
			<<_triangles[i].v3<<std::endl;
	}
	file<<std::endl;

	file<<"CELL_TYPES "<<_triangles.size()<<std::endl;
	for(int i=0; i<_triangles.size(); ++i){
		file<<"5"<<std::endl;
	}
	file<<std::endl;

	file.close();

	return 1;	
}
