// Motor test and print calibration data to serial monitor

#include <Balboa32U4.h>
#include <Wire.h>
#include <LSM6.h>

LSM6 imu;
Balboa32U4Motors motors;
Balboa32U4Encoders encoders;

const uint8_t UPDATE_TIME_MS = 10;
const uint8_t CALIBRATION_ITERATIONS = 100;

uint32_t startTime_ms;
uint32_t prevTime_ms;
uint32_t curTime_ms;
float curTime_s;

int32_t gYZero;
int32_t angleIntegrated;

float angleGravity;
float angleGravityTemp;

int32_t gY;
int32_t aX;
int32_t aZ;

float angle_rad;
float alpha = .8;

void setup()
{
  Serial.begin(115200);
  gyroSetup();

  startTime_ms = millis();
  prevTime_ms = startTime_ms;

  Serial.println("time_s angle_rad angle_deg gravity_deg");
  
}

void loop()
{
  curTime_ms = millis();
  curTime_s = ((float) curTime_ms) / 1000.0;

  if (curTime_ms - prevTime_ms >= UPDATE_TIME_MS)
  {
    imu.read();
    integrateGyro();

    float angle_deg = angle_rad * 180.0 / 3.14159;
    float gravity_deg = angleGravity * 180.0 / 3.14159;

    Serial.print(curTime_s,3);
    Serial.print(" ");

    Serial.print(angle_rad,6);
    Serial.print(" ");

    Serial.print(angle_deg,3);
    Serial.print(" ");

    Serial.println(gravity_deg,3);

    Serial.print("aX = "); Serial.print(imu.a.x);
    Serial.print("  aY = "); Serial.print(imu.a.y);
    Serial.print("  aZ = "); Serial.println(imu.a.z);


    prevTime_ms = curTime_ms;
  }
}

void gyroSetup()
{
  Wire.begin();
  if (!imu.init())
  {
    while(true)
    {
      Serial.println("Failed to detect and initialize IMU!");
      delay(200);
    }
  }

  imu.enableDefault();
  imu.writeReg(LSM6::CTRL2_G,  0b01011000);
  imu.writeReg(LSM6::CTRL1_XL, 0b01011000);

  delay(1000);

  int32_t totalY = 0;
  for (int i = 0; i < CALIBRATION_ITERATIONS; i++)
  {
    imu.read();
    totalY += imu.g.y;
    delay(1);
  }

  gYZero = totalY / CALIBRATION_ITERATIONS;

  angleIntegrated = (int32_t)(atan2f((double)imu.a.z,(double)imu.a.x)*57295);
}

void integrateGyro()
{
  gY = (imu.g.y - gYZero) / 29;

  aX = imu.a.x;
  aZ = imu.a.z;

  angleIntegrated += gY * (curTime_ms - prevTime_ms);

  angle_rad = ((float)angleIntegrated)/1000/180*3.14159;

  angleGravity = atan2f((double)aZ,(double)aX);

  angleGravityTemp = angleGravity;

  while(angleGravityTemp-angle_rad >= 3.14159)
    angleGravityTemp -= 6.28319;

  while(angleGravityTemp-angle_rad <= -3.14159)
    angleGravityTemp += 6.28319;

  angle_rad = (1-alpha)*angle_rad + alpha*angleGravityTemp;
}