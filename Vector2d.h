///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: CVector2d                                                     //
//                                                                       //
//  2D Vector to represent points or directions.  Each component of the  //
//  vector is stored as a doubleing point number.                        //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _VECTOR_2D_
#define _VECTOR_2D_

class CVector2d
{

private :

protected:

 double vec[2];  // Storage for the vector components

public :

  // Constructors
  CVector2d() { }
  CVector2d(const double x, const double y)
    { vec[0] = x; vec[1] = y; }
   CVector2d(const double v[2])
    { vec[0] = v[0]; vec[1] = v[1];}
  CVector2d(const CVector2d &vector)  { Set(vector); }
  CVector2d(const CVector2d *pVector) { Set(pVector); }
  CVector2d(const CVector2d &a, const CVector2d& b)
    { Set(b - a); }
  CVector2d(const CVector2d *a, const CVector2d *b)
    { Set(*b - *a); }
  
  virtual ~CVector2d() { }

  // Debug
  void Trace() const;

  // Data setting
  void Clear()
    { vec[0] = 0.f; vec[1] = 0.f;  }
  void Set(const CVector2d *pVector) { Set(pVector->GetArray()); }
  void Set(const CVector2d& vector)  { Set(vector.GetArray()); }
  void Set(const double x, const double y)
    { vec[0] = x; vec[1] = y;  }
  void Set(const double v[1])
    { vec[0] = v[0]; vec[1] = v[1];}
  void Set(const CVector2d& a, const CVector2d& b)
    { Set(b - a); }
  void Set(const CVector2d *a, const CVector2d *b)
    { Set(*b - *a); }

  // Data Access
  const double* GetArray() const { return vec; }
  void         Get(double& x, double& y) const;

  // Per coordinate (explicit inline functions)
  void x(double newX) { vec[0] = newX; }
  void y(double newY) { vec[1] = newY; }

  // Data access (explicit inline functions)
  double x() const { return (vec[0]); }
  double y() const { return (vec[1]); }

  // Data access using indices
  double&       operator[](int i)       { return (vec[i]); }
  const double& operator[](int i) const { return (vec[i]); }

  // Operators
  CVector2d& operator+=(const CVector2d& rVector);
  CVector2d& operator+=(const CVector2d* pVector);
  CVector2d& operator-=(const CVector2d& rVector);
  CVector2d& operator-=(const CVector2d* pVector);
  CVector2d& operator*=(double d);
  CVector2d& operator/=(double d)
    { return *this *= (1.f/d); }

  // Nondestructive unary negation - returns a new vector
  CVector2d  operator -() const;

  // Binary operators
  friend CVector2d operator+(const CVector2d& u, const CVector2d& v);
  friend CVector2d operator-(const CVector2d& u, const CVector2d& v);
  friend CVector2d operator*(double s,            const CVector2d& u);
  friend CVector2d operator*(const CVector2d& u, double s)
    { return s * u; }
  friend CVector2d operator/(const CVector2d& u, double s)
    { return (1.f/s) * u; }
  friend double operator^(const CVector2d& u, const CVector2d& v);
  friend int       operator==(const CVector2d& v1, const CVector2d& v2);
  friend int       operator!=(const CVector2d& v1, const CVector2d& v2)
    { return !(v1 == v2); }

  int Equals(const CVector2d& v, double tolerence) const;

  inline double     Dot(const CVector2d& v) const;
  inline double     Dot(const CVector2d* pV) const;
  double        Cross(const CVector2d& v) const;
  double        Cross(const CVector2d* pV) const;

  // Misc
  double Self_Operator(CVector2d *vector);
  double Normalize();
  double Normalize(double value);
  double Length() const;
  double LengthSquared() const;
  int IsCollinear(CVector2d *pVector) const;
  int IsCollinear(CVector2d &vector) const;
	void Negate();
	//CVector2d Rotate(double angle,CVector2d Around); 
	//CVector2d Projection(const CVector2d* pV) const;
	//double Elem_sum();
};


//********************************************
// Dot
//********************************************
double
CVector2d::Dot(const CVector2d& v) const
{
  return (x() * v.x() +
	  y() * v.y());
}

//********************************************
// Dot
//********************************************
double
CVector2d::Dot(const CVector2d* pV) const
{
  return (x() * pV->x() +
	  y() * pV->y());
}

#endif // _VECTOR_3D_
