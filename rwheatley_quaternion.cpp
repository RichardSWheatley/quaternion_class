/********************************************************
 ********************************************************
 * Quaternion class
 * By Richard Wheatley (richardswheatley@gmail.com)
 * Original: 11/23/2016  
 * Revised:  
 *   
 * Dependancies: 
 * 
 * � 2016, This code is provided "as is" with no warranty
 *   You may use it freely as long as credit is given to
 *   Richard S Wheatley in the application in which it is
 *   used
 ********************************************************
 ********************************************************/

 #include <math.h>
// Constants
const double PI = 3.141592653589793;
const double halfPI = 1.570796326794897;
const double twoPI = 6.283185307179586;

class vector3_t{
private:
	double vx, vy, vz;

protected:

	
public:
	// a simple destructor.
	virtual ~vector3_t();
	
	// Constructors with and without initialized values
	vector3_t(void){
		vx = 0.0; vy = 0.0; vz = 0.0;
	}

	vector3_t(double a, double b, double c){
		vx = a; vy = b; vz = c;
	}

	double getVX(){ return(vx);}
	double getVY(){ return(vx);}
	double getVZ(){ return(vx);}
};




class quat_t {
private:

	/********************************************************
	 ********************************************************
	 * The 'w' below is the "real" part of the quaternion
	 * More information about quaternions can be found by
	 * serching for "Hamliton's Quaternions"
	 ********************************************************
	 ********************************************************/
	double qw;

	/********************************************************
	 ********************************************************
	 * The 'qx', 'qy', 'qz' below are the "imaginary" parts of
	 * the quaternion. More information about quaternions can
	 * be found by serching for "Hamliton's Quaternions"
	 ********************************************************
	 ********************************************************/
	double qx, qy, qz;

protected:
	// Decontructors for this quaternion class
	virtual ~quat_t();

public:
	// Contructors for this quaternion class
	quat_t() {}
	quat_t(double real, double i, double j, double k) { qw = real; qx = i; qy = j; qz = k; }
	quat_t(double real, vector3_t values) { qw = real; qx = values.getVX(); qy = values.getVY(); qz = values.getVZ(); }
	quat_t(const quat_t& q){ qw = q.qw; qx = q.qx; qy = q.qy; qz = q.qz; }

	void setQW(double _qw){ qw = _qw; }
	void setQX(double _qx){ qx = _qx; }
	void setQY(double _qy){ qy = _qy; }
	void setQZ(double _qz){ qz = _qz; }

	double getQW(){ return qw; }
	double getQX(){ return qx; }
	double getQY(){ return qy; }
	double getQZ(){ return qz; }

	// basic operations
	quat_t &operator =(const quat_t &q)	{
		qw = q.qw;
		qx = q.qx;
		qy = q.qy;
		qz = q.qz;
		return *this;
	}

	const quat_t operator +(const quat_t &q) const {
		return quat_t(qw + q.qw, qx + q.qx, qy + q.qy, qz + q.qz);
	}

	const quat_t operator -(const quat_t &q) const {
		return quat_t(qw - q.qw, qx - q.qx, qy - q.qy, qz - q.qz);
	}

	const quat_t operator *(const quat_t &q) const {
		return quat_t(qw*q.qw - qx*q.qx - qy*q.qy - qz*q.qz,
						  qx*q.qw + qw*q.qx - qz*q.qy + qy*q.qz,
						  qy*q.qw + qz*q.qx + qw*q.qy - qx*q.qz,
						  qz*q.qw - qy*q.qx + qx*q.qy + qw*q.qz);
	}

	const quat_t operator /(const quat_t &q) const {	
			quat_t Q(q); 
			Q.inverse(); 
			return *this * Q;
	}

	const quat_t operator *(double scalar_value) const {
		return quat_t(qw*scalar_value, qx*scalar_value, qy*scalar_value, qx*scalar_value);
	}

	const quat_t operator /(double scalar_value) const {
		return quat_t(qw / scalar_value, qx / scalar_value, qy / scalar_value, qz / scalar_value);
	}

	const quat_t operator -() const	{
			return quat_t(-qw, -qx, -qy, -qz);
	}
	
	const quat_t &operator +=(const quat_t &q) {
		qw += q.qw;
		qx += q.qx;
		qy += q.qy;
		qz += q.qz; 
		return *this;
	}

	const quat_t &operator -=(const quat_t &q) {
		qw -= q.qw;
		qx -= q.qx;
		qy -= q.qy;
		qz -= q.qz;
		return *this;
	}

	const quat_t &operator *=(const quat_t &q) {	
		double w = qw, x = qx, y = qy, z = qz, sn = qw*q.qw - qx*q.qx - qy*q.qy - qz*q.qz;
		qw = sn;
		qx = w*q.qx + x*q.qw + y*q.qz - z*q.qy;
		qy = w*q.qy - x*q.qz + y*q.qw + z*q.qx;
		qz = w*q.qz + x*q.qy - y*q.qx + z*q.qw;
		return *this;
	}
	
	const quat_t &operator *=(double scalar_value) {
		qw *= scalar_value;
		qx *= scalar_value;
		qy *= scalar_value;
		qz *= scalar_value;
		return *this;
	}

	const quat_t &operator /=(double scalar_value) {
		qw /= scalar_value;
		qx /= scalar_value;
		qy /= scalar_value;
		qz /= scalar_value;
		return *this;
	}
	
	// gets the length of this quat_t
	double get_length() const {
		return sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
	}

	// gets the squared length of this quat_t
	double get_length_squared() const {
		return (qw*qw + qx*qx + qy*qy + qz*qz);
	}

	// normalizes this quat_t
	void normalize() {
		*this /= get_length();
	}

	// returns the normalized version of this quat_t
	quat_t normalized() const {
		return *this / get_length();
	}

	// Forces this quat_t to be its conjugate
	void conjugate() {
		qx = -qx; qy = -qy; qz = -qz;
	}

	// Forces this quat_t to be its inverse
	void inverse() {
		conjugate();
		*this /= get_length_squared();
	}
	
	// computes the dot product of 2 quat_ts
	static inline double dot(const quat_t &q_left, const quat_t &q_right) {
		return (q_left.qw*q_right.qw + q_left.qx*q_right.qx + q_left.qy*q_right.qy + q_left.qz*q_right.qz);
	}
};


// linear quat_t interpolation
static quat_t lerp(const quat_t &q_left, const quat_t &q_right, double time) {

}

// spherical linear interpolation
static quat_t slerp(const quat_t &q_left, const quat_t &q_right, double time)
{

}

// Elevation Quaternion from Normalized Accel data
//   Accel data must be normalized or the
//   equation asin(-ax) won't work.
void pitch_quaternion(vector3_t &a, quat_t &quat_ptr){
	double angle_pitch = atan2(-a.getVX(), sqrt(a.getVY()*a.getVY() + a.getVZ()*a.getVZ()));

	quat_ptr.setQW(cos(angle_pitch/2.0));
	quat_ptr.setQX(0.0);
	quat_ptr.setQY(sin(angle_pitch/2.0));
	quat_ptr.setQZ(0.0);
}

// Roll Quaternion from Normalized Accel data
//   Accel data must be normalized or the
//   equation angle equations won't work.
void roll_quaternion(vector3_t &a, quat_t &quat_ptr){
	double angle_roll = atan2(a.getVY(), a.getVZ());

	quat_ptr.setQW(cos(angle_roll/2.0));
	quat_ptr.setQX(sin(angle_roll/2.0));
	quat_ptr.setQY(0.0);
	quat_ptr.setQZ(0.0);
}

// Azimuth Quaternion from Normalized mag data
void azimuth_quaternion(quat_t &q_e, quat_t &q_r, quat_t &q_m, quat_t &quat_ptr){
	quat_ptr = (quat_t)(q_e * q_r * q_m / q_r / q_e);
	double squared_magnitude = sqrt(quat_ptr.getQX() * quat_ptr.getQX() + quat_ptr.getQY() * quat_ptr.getQY());
	double angle_az = atan2(quat_ptr.getQY()/squared_magnitude, quat_ptr.getQX()/squared_magnitude);

	quat_ptr.setQW(cos(angle_az/2.0));
	quat_ptr.setQX(0.0);
	quat_ptr.setQY(0.0);
	quat_ptr.setQZ(sin(angle_az/2.0));
}

// Total Quaternion from Normalized mag data
void total_rotation_quaternion(quat_t &q_a, quat_t &q_e, quat_t &q_r, quat_t &quat_ptr){
	quat_ptr = q_a * q_e * q_r;
}

void quat_vect_multiply(quat_t &q, vector3_t &v, quat_t &quat_ptr){
	quat_ptr.setQW(0.5*(-q.getQX() * v.getVX() - q.getQY() * v.getVY() - q.getQZ() * v.getVZ()));
	quat_ptr.setQX(0.5*( q.getQW() * v.getVX() - q.getQZ() * v.getVY() + q.getQY() * v.getVZ()));
	quat_ptr.setQY(0.5*( q.getQZ() * v.getVX() + q.getQW() * v.getVY() - q.getQX() * v.getVZ()));
	quat_ptr.setQZ(0.5*(-q.getQY() * v.getVX() + q.getQX() * v.getVY() + q.getQW() * v.getVZ()));
}

"# quaternion_class" 
