// 10 bit ADC
//A code to read a voltage 0-5V as [0,1023]
//Use a potentiometer's ends to 5V and 0V
//Connect the middle pot pin to A2 on the Arduino. 



void setup() {

pinMode(A2, INPUT);

Serial.begin(9600);
}

void loop() {

int x = analogRead(A2);

Serial.print("Value is ");
Serial.println(x);

//can be read from Serial Monitor
//Magnifying glass in upper right corner
//Match Baudrate (9600)

delay(50);
}

//Now let's use a photoresistor!
//Did Drew remember the 1k Resistors? 
