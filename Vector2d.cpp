//********************************************
// Vector2d.cpp
//********************************************
// class CVector2d
//********************************************
// mmeyer@gg.caltech.edu
// Created :  09/07/00
// Modified : 09/07/00
// Modified by beibei for 2D vector
//********************************************


#include "Vector2d.h"

#include <math.h>
#include <stdio.h>


//////////////////////////////////////////////
// CONSTRUCTORS
//////////////////////////////////////////////


//////////////////////////////////////////////
// DATA
//////////////////////////////////////////////

//********************************************
// Get
//  Get the 2 components of the vector
//********************************************
void
CVector2d::Get(double& x, double& y) const
{
  x = vec[0];
  y = vec[1];
}

//********************************************
// Trace
//********************************************
void
CVector2d::Trace() const
{/*
  TRACE("\n");
  TRACE("** Vector **\n");
  TRACE("Address      : %x\n",this);
  TRACE("Coordinates : (%g %g %g)\n",vec[0],vec[1],vec[2]);*/
}

//////////////////////////////////////////////
// OPERATORS
//////////////////////////////////////////////

//********************************************
// Operator +=
//********************************************
CVector2d&
CVector2d::operator+=(const CVector2d& rVector)
{
  vec[0] += rVector.x();
  vec[1] += rVector.y();
  return *this;
}

//********************************************
// Operator +=
//********************************************
CVector2d&
CVector2d::operator+=(const CVector2d* pVector)
{
  vec[0] += pVector->x();
  vec[1] += pVector->y();
  return *this;
}

//********************************************
// Operator -=
//********************************************
CVector2d&
CVector2d::operator-=(const CVector2d& rVector)
{
  vec[0] -= rVector.x();
  vec[1] -= rVector.y();
  return *this;
}

//********************************************
// Operator -=
//********************************************
CVector2d&
CVector2d::operator-=(const CVector2d* pVector)
{
  vec[0] -= pVector->x();
  vec[1] -= pVector->y();
  return *this;
}

//********************************************
// Operator *=
//********************************************
CVector2d&
CVector2d::operator*=(double d)
{
  vec[0] *= d;
  vec[1] *= d;
  return *this;
}

//********************************************
// Operator -
//  Nondestructive unary -
//  Returns a new vector.
//********************************************
CVector2d
CVector2d::operator -() const
{
  return CVector2d(-vec[0],-vec[1]);
}

//********************************************
// Operator + 
//********************************************
CVector2d
operator+(const CVector2d& u, const CVector2d& v)
{
  return CVector2d(u.vec[0]+v.vec[0],u.vec[1]+v.vec[1]);
}

//********************************************
// Operator -
//********************************************
CVector2d
operator-(const CVector2d& u, const CVector2d& v)
{
  return CVector2d(u.vec[0]-v.vec[0],u.vec[1]-v.vec[1]);
}

//********************************************
// Operator * 
//********************************************
CVector2d
operator*(double s, const CVector2d& u)
{
  return CVector2d(u.vec[0] * s, u.vec[1] * s);
}

//********************************************
// Operator ^
//  Returns the cross product of u and v.
//********************************************
double
operator^(const CVector2d& u, const CVector2d& v)
{
  return (u.vec[0] * v.vec[1] - u.vec[1] * v.vec[0]);
}

//********************************************
// Operator ==
//********************************************
int
operator==(const CVector2d& v1, const CVector2d& v2)
{
  return (v1.vec[0] == v2.vec[0] &&
	  v1.vec[1] == v2.vec[1]);
}

//********************************************
// Equals
//  Determines if two vectors are equal
//  within a tolerence (squared distance).
//********************************************
int
CVector2d::Equals(const CVector2d& v, double tolerence) const
{
  CVector2d diff = *this - v;

  return diff.LengthSquared() <= tolerence;
}


//////////////////////////////////////////////
//////////////////////////////////////////////
// PROCESSING
//////////////////////////////////////////////
//////////////////////////////////////////////


//********************************************
// Cross
//********************************************
double
CVector2d::Cross(const CVector2d& v) const
{
  return (x() * v.y() - y() * v.x());
}

//********************************************
// Cross
//********************************************
double
CVector2d::Cross(const CVector2d* pV) const
{
  return  (x() * pV->y() - y() * pV->x());
}


//********************************************
// Normalize
//********************************************
double
CVector2d::Normalize()
{
  double len = Length();
  if(len != 0.0f)
    (*this) *= (1.0/len);
  else
    Set(0.f,0.f);

  return len;
}

//********************************************
// Normalize
//********************************************
double
CVector2d::Normalize(double value)
{
  double len = Length();
  if(len != 0.0f)
    (*this) *= (value/len);
  else
    Set(0.f,0.f);

  return len;
}

//********************************************
// LengthSquared
//********************************************
double
CVector2d::LengthSquared() const
{
  return ( (double)vec[0]*(double)vec[0]
	 + (double)vec[1]*(double)vec[1]);
}
	
//********************************************
// Length
//********************************************
double
CVector2d::Length() const
{
  return sqrt( (double)vec[0]*(double)vec[0]
	     + (double)vec[1]*(double)vec[1]);
}
	
//********************************************
// IsCollinear
//********************************************
int
CVector2d::IsCollinear(CVector2d *pVector) const
{
  double x = pVector->x() / vec[0];
  double y = pVector->y() / vec[1];
  return ((x == y));
}

//********************************************
// IsCollinear
//********************************************
int
CVector2d::IsCollinear(CVector2d &vector) const
{
  double x = vector.x() / vec[0];
  double y = vector.y() / vec[1];

  return ((x == y));
}

//********************************************
// Negate
//********************************************
void
CVector2d::Negate()
{
  vec[0] = -vec[0];
  vec[1] = -vec[1];
 
}

////********************************************
//// Rotate this vector around pAround by angle
//// by Haeyoung Lee
////********************************************
//CVector2d CVector2d::Rotate(double angle, 
//														CVector2d Around) 
//{
//	double f1, f2, f3;
//	CVector2d t1, t2;
//	
//	
//	f1 = (double)cos((double)angle);
//	f2 = (double)sin((double)angle);
//	t1 = Projection(&Around);
//	t2 = Around.Cross(this);
//	f3 = Dot(Around);
//	
//	return CVector2d((double)(f1*t1.x()+f2*t2.x()+f3*Around.x()),
//		(double)(f1*t1.y()+f2*t2.y()+f3*Around.y()),
//		(double)(f1*t1.z()+f2*t2.z()+f3*Around.z()));
//	
//}
//
////********************************************
//// Projection
//// by Haeyoung Lee
////********************************************
//CVector2d    
//CVector2d::Projection(const CVector2d* pV) const
//{
//  double alpha = Dot(pV)/pV->Dot(pV);
//	return CVector2d(x()-alpha* pV->x(), 
//		               y()-alpha*pV->y(),
//		               z()-alpha*pV->z());
//}
//
//
//double 
//CVector2d::Self_Operator(CVector2d *vector){
//	return(x()*vector->y()-y()*vector->x());
//}
//
//double CVector2d::Elem_sum(){
//	return (x() - y() + z());
//
//}


// ** EOF **
