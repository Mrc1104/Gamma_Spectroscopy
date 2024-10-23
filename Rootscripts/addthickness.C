#include "include/GPHYS_Quantity.h"

using GQ = GPHYS_Quantity;
void addthickness()
{
	GQ a0(0,0);
	GQ a1(1.29,0.01);
	GQ a2(1.20, 0.01);
	GQ a3(1.28, 0.01);

	cout << "Aluminum:" << endl;
	cout << "-------------------------------------" << endl;
	cout << "0 Plate: " << a0 << endl;
	cout << "1 Plate: " <<  a0 + a1 << endl;
	cout << "2 Plate: " << a0 + a1 + a2 << endl;
	cout << "3 Plate: " << a0 + a1 + a2 + a3 << endl;

	GQ c1(1.25, 0.01);
	GQ c2(1.29, 0.01);
	GQ c3(1.29, 0.01);

	cout << "Copper:" << endl;
	cout << "-------------------------------------" << endl;
	cout << "0 Plate: " << a0 << endl;
	cout << "1 Plate: " <<  a0 + c1 << endl;
	cout << "2 Plate: " <<  a0 + c1 + c2 << endl;
	cout << "3 Plate: " <<  a0 + c1 + c2 + c3 << endl;


}
