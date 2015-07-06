//Code via James @ MakeUseOf.com

int data = 11;
int clock = 12;
int latch = 8;
int x = 0;

void setup()
{
   // set the three control pins to output
  pinMode(data, OUTPUT);
  pinMode(clock, OUTPUT);  
  pinMode(latch, OUTPUT);  
    Serial.begin(9600); // so we can send debug messages to serial monitor
}

void loop(){
   
    outputBytes(); // our basic output which writes 8-bits to show how a shift register works. 
    //outputIntegers(); // sends an integer value as data instead of bytes, effectively counting in binary. 
}
 
void outputIntegers(){
     for (int i=0;i<256;i++){
        digitalWrite(latch, LOW);     
        Serial.println(i);  // Debug, sending output to the serial monitor
        shiftOut(data, clock, MSBFIRST, i); 
        digitalWrite(latch, HIGH);   
        delay(100);    
     } 
}

void outputBytes(){
    /* Bytes, or 8-bits, are represented by a B followed by 8 0 or 1s. 
        In this instance, consider this to be like an array that we'll use to control
        the 8 LEDs. Here I've started the byte value as 00000001
    */    
        
    byte dataValues = B00000001; // change this to adjust the starting pattern
    
    /* In the for loop, we begin by pulling the latch low, 
        using the shiftOut Arduino function to talk to the shift register, 
        sending it our byte of dataValues representing the state of the LEDs
        then pull the latch high to lock those into place.
        
        Finally, we shift the bits one place to the left, meaning the next iteration
        will turn on the next LED in the series.
        
        To see the exact binary value being sent, check the serial monitor.
    */
    
    for (int i=0;i<8;i++){
      digitalWrite(latch, LOW);     
      Serial.println(dataValues, BIN);  // Debug, sending output to the serial monitor
      shiftOut(data, clock, MSBFIRST, dataValues); 
      digitalWrite(latch, HIGH);   
      dataValues = dataValues << 1; // Shift the bits one place to the left -  change to >> to adjust direction
      delay(100);    
    }
 
}
