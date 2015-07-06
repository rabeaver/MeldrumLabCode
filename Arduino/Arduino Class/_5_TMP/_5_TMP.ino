//Using TMP36 3-pin thermometer
//Ground = right side (as reading front)
//Power = left side
//Data out is middle pin, wire to A0


void setup() {
  
pinMode(A0, INPUT);
Serial.begin(9600);
}

void loop() {
int reading = analogRead(A0); //read V(T)

/*Serial.print("V is: ");
Serial.println(reading);
*/

float vtemp = reading *5.0; 
vtemp/=1024.0; //convert to a 5v scale
//Serial.println(vtemp);

float TempC = (vtemp-.5)*100; //10mV/deg, offset .5v
float TempF = TempC*(9./5.) + 32.;

//Serial.print("Temp in C is: ");
Serial.print(TempC);
Serial.println();

//Serial.print("Temp in F is: ");
//Serial.println(TempF);

delay(60000);
}

//Later, output this to the LCD Screen!
