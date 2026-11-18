#include <Wire.h>
#include "SparkFun_BMA400_Arduino_Library.h"

#define SDA_PIN 24
#define SCL_PIN 25
#define ADDR    0x15

#define M       5
#define FMIN    0.8
#define EPS     3.0
#define NU      1.0
#define MU      1.0
#define EPSH    2.0
#define FS      200.0
#define SETTLE  2.0

BMA400 accel;

float phi[M], alpha[M], omega, beta;
float phiF[M], alphaF[M], omegaF, betaF;
float hx, hy, hw;
float dcLevel, ampLevel;
bool started;
unsigned long nextTick;
unsigned long n;

const float dt   = 1.0 / FS;
const float wMin = 2.0 * PI * 0.2;
const float wMax = 2.0 * PI * 5.0;

void resetState() {
  for (int i = 0; i < M; i++) {
    phi[i] = 0; alpha[i] = 0;
    phiF[i] = 0; alphaF[i] = 0;
  }
  omega = 2.0 * PI * FMIN;
  omegaF = omega;
  beta = 0; betaF = 0;
  hx = 1.0; hy = 0.0; hw = omega;
  ampLevel = 1.0;
  started = false;
  n = 0;
}

float update(float *ph, float *al, float *w, float *bt, float target, bool learn, bool adapt) {
  float est = *bt;
  for (int i = 0; i < M; i++) est += al[i] * sin(ph[i]);

  float err = target - est;

  float dPhi[M], dAl[M];
  for (int i = 0; i < M; i++) {
    dPhi[i] = (i + 1) * (*w);
    dAl[i]  = 0;
    if (learn) {
      dPhi[i] += EPS * err * cos(ph[i]);
      dAl[i]   = NU * err * sin(ph[i]);
    }
  }

  for (int i = 0; i < M; i++) {
    ph[i] += dPhi[i] * dt;
    if (ph[i] >= 2.0 * PI) ph[i] -= 2.0 * PI;
    al[i] += dAl[i] * dt;
  }

  if (learn) {
    if (adapt) {
      *w += EPS * err * cos(ph[0]) * dt;
      if (*w < wMin) *w = wMin;
      if (*w > wMax) *w = wMax;
    }
    *bt += NU * err * dt;
  }

  return est;
}

void updateHopf(float f, bool learn) {
  float r = sqrt(hx * hx + hy * hy);
  if (r < 1e-6) r = 1e-6;

  float drive = learn ? f : 0.0;

  float dx = (MU - r * r) * hx - hw * hy + EPSH * drive;
  float dy = (MU - r * r) * hy + hw * hx;
  float dw = -EPSH * drive * hy / r;

  hx += dx * dt;
  hy += dy * dt;
  hw += dw * dt;
  if (hw < wMin) hw = wMin;
  if (hw > wMax) hw = wMax;
}

void printHeader() {
  Serial.println("t,theta_d,theta_hat,theta_fixed,freq,alpha1,phi1,x,y,freq_hopf");
}

void setup() {
  Serial.begin(115200);
  delay(1000);

  Wire.setSDA(SDA_PIN);
  Wire.setSCL(SCL_PIN);
  Wire.begin();
  Wire.setClock(400000);

  while (accel.beginI2C(ADDR) != BMA400_OK) {
    Serial.println("BMA400 not found");
    delay(1000);
  }
  accel.setRange(BMA400_RANGE_2G);
  accel.setODR(BMA400_ODR_400HZ);

  resetState();
  printHeader();
  nextTick = micros();
}

void loop() {
  while (Serial.available()) {
    if (Serial.read() == 'r') {
      resetState();
      printHeader();
    }
  }

  if ((long)(micros() - nextTick) < 0) return;
  nextTick += (unsigned long)(1000000.0 / FS);

  if (accel.getSensorData() != BMA400_OK) return;

  float gx = accel.data.accelX;
  float gy = accel.data.accelY;
  float gz = accel.data.accelZ;
  float mag = sqrt(gx * gx + gy * gy + gz * gz);

  if (!started) { dcLevel = mag; started = true; }
  dcLevel += 2.0 * PI * 0.15 * dt * (mag - dcLevel);

  float theta = mag - dcLevel;
  ampLevel += 2.0 * PI * 0.20 * dt * (fabs(theta) - ampLevel);
  if (ampLevel < 0.02) ampLevel = 0.02;
  theta /= ampLevel;

  n++;
  bool learn = (n > (unsigned long)(SETTLE * FS));

  float est  = update(phi,  alpha,  &omega,  &beta,  theta, learn, true);
  float estF = update(phiF, alphaF, &omegaF, &betaF, theta, learn, false);
  updateHopf(theta, learn);

  Serial.print(n * dt, 4);            Serial.print(',');
  Serial.print(theta, 4);             Serial.print(',');
  Serial.print(est, 4);               Serial.print(',');
  Serial.print(estF, 4);              Serial.print(',');
  Serial.print(omega / (2.0 * PI), 4); Serial.print(',');
  Serial.print(alpha[0], 4);          Serial.print(',');
  Serial.print(phi[0], 4);            Serial.print(',');
  Serial.print(hx, 4);                Serial.print(',');
  Serial.print(hy, 4);                Serial.print(',');
  Serial.println(hw / (2.0 * PI), 4);
}
