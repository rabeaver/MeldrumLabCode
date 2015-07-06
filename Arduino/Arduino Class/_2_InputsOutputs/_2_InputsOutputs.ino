// 10 bit ADC
//A code to read a voltage 0-5V as [0,1023]
//Use a potentiometer's ends to 5V and 0V
//Connect the middle pot pin to A2 on the Arduino. 
//Connect an LED to pin 3


void setup() {
//setup code here, to run once:
pinMode(A2, INPUT);
pinMode(3, OUTPUT);
Serial.begin(9600);
}

void loop() {
// main code here, to run repeatedly:
int minV;
int maxV;
int x = analogRead(A2);
int y = map(x, 0, 1023, 0, 255);
analogWrite(3, y);
Serial.print("Value is ");
Serial.println(x);

//can be read from Serial Monitor
//Magnifying glass in upper right corner
//Match Baudrate

delay(50);
}

//Now let's use a photoresistor!
