#include <Python.h>

double *
_RandomAngle(double Jx0,double Jy0,double Jz0,double angle,double random)
{
	static double r[3];
	double modulo = sqrt(Jx0*Jx0+Jy0*Jy0+Jz0*Jz0);
	if (modulo < 1e-10){
		r[0] = 0;
		r[1] = 0;
		r[2] = 0;
		return  r;
	}
	double Jx,Jy,Jz,angulo_t,angulo_phi,Jx2,Jy2,Jz2;
	Jz0 = Jz0/modulo;
	Jy0 = Jy0/modulo;
	Jx0 = Jx0/modulo;
	angulo_t = acos(Jz0)+angle ;
	angulo_phi = atan2(Jy0,Jx0);
	Jz = cos(angulo_t);
	Jy = sin(angulo_t)*sin(angulo_phi);
	Jx = sin(angulo_t)*cos(angulo_phi);
	double CosA = cos(random);
	double SinA = sin(random);
	Jx2 = (CosA+Jx0*Jx0*(1-CosA))*Jx       +   (Jx0*Jy0*(1-CosA)-Jz0*SinA)*Jy      + (Jx0*Jz0*(1-CosA)+Jy0*SinA)*Jz;
	Jy2 = (Jx0*Jy0*(1-CosA)+Jz0*SinA)*Jx   +     (CosA+Jy0*Jy0*(1-CosA))*Jy        + (Jy0*Jz0*(1-CosA)-Jx0*SinA)*Jz;
	Jz2 =  (Jx0*Jz0*(1-CosA)-Jy0*SinA)*Jx  +    (Jz0*Jy0*(1-CosA)+Jx0*SinA)*Jy     + (CosA+Jz0*Jz0*(1-CosA))*Jz;

	
	r[0] = Jx2*modulo;
	r[1] = Jy2*modulo;
	r[2] = Jz2*modulo;
	return  r;
}
static PyObject*
RandomAngle(PyObject* self, PyObject* args)
{
    double Jx0;
    double Jy0;
    double Jz0;
    double angle;
    double random;


    if (!PyArg_ParseTuple(args, "ddddd", &Jx0,&Jy0,&Jz0,&angle,&random))
        return NULL;

    return Py_BuildValue("ddd", _RandomAngle(Jx0,Jy0,Jz0,angle,random)[0],_RandomAngle(Jx0,Jy0,Jz0,angle,random)[1],_RandomAngle(Jx0,Jy0,Jz0,angle,random)[2]);
}
static PyMethodDef Prop[] = {
    {"RandomAngle", RandomAngle, METH_VARARGS, "Random Or. a 3d vector in an angle alpha"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initRandomAngle(void)
{
    (void) Py_InitModule("RandomAngle", Prop);
}