// Timestamped Temp and RH monitor
// From ladyada code

#include "SD.h"
#include "Wire.h"
#include "RTClib.h"
#include "DHT.h"
#include "SPI.h"


#define DHTPIN 2     // what pin we're connected to for T/RH sensor
#define LOGDELAY 900000 // time delay for logging
#define ECHO_TO_SERIAL 0 //output to display as well
#define DHTTYPE DHT22   // define RH/Temp sensor type DHT 22  (AM2302)
#define redLEDpin 3
#define greenLEDpin 4

DHT dht(DHTPIN, DHTTYPE); // Initialize DHT sensor for normal 16mhz Arduino
RTC_DS1307 RTC; // define the Real Time Clock object
const int chipSelect = 10; // for the data logging shield, we use digital pin 10 for the SD cs line

// the logging file
File logfile;



  
void setup(void) {
  Serial.begin(9600); 

  
  // initialize the SD card
  //Serial.print("Initializing SD card...");
  
//  // make sure that the default chip select pin is set to
//  // output, even if you don't use it:
//  pinMode(chipSelect, OUTPUT);
  
  // see if the card is present and can be initialized:
  if (!SD.begin(chipSelect)) {
    Serial.println("Card failed, or not present");
    // don't do anything more:
    return;
  }
  //Serial.println("card initialized.");
  

  // create a new file
  char filename[] = "LOGGER00.txt";
  for (uint8_t i = 0; i < 100; i++) {
    filename[6] = i/10 + '0';
    filename[7] = i%10 + '0';
    if (! SD.exists(filename)) {
      // only open a new file if it doesn't exist
      logfile = SD.open(filename, FILE_WRITE); 
      break;  // leave the loop!
    }
  }
  
 // open the file. note that only one file can be open at a time,
  // so you have to close this one before opening another.
//  logfile = SD.open(filename, FILE_WRITE);
  
  
//  if (logfile) {
//    Serial.print("Writing to test.txt...");
//    logfile.println("testing 1, 2, 3.");
//    // close the file:
//    logfile.close();
//    Serial.println("done.");
//  } else {
//    // if the file didn't open, print an error:
//    Serial.println("error opening test.txt");
//  }
  
  
  //Serial.print("Logging to: ");
 // Serial.println(filename);
  
    Wire.begin();  
  if (!RTC.begin()) {
    logfile = SD.open(filename, FILE_WRITE);
    logfile.println("RTC failed");
    logfile.close();
#if ECHO_TO_SERIAL
    Serial.println("RTC failed");
#endif  //ECHO_TO_SERIAL
  }
  
  logfile = SD.open(filename, FILE_WRITE);
  logfile.println("date,time,rh,temp_C");  
  logfile.close();  
#if ECHO_TO_SERIAL
  Serial.println("date,time,rh,temp_C");
  
#endif

  pinMode(redLEDpin, OUTPUT);
  pinMode(greenLEDpin, OUTPUT);
}


void loop(void) {
  // Wait a few seconds between measurements.

  delay(LOGDELAY);
 logfile = SD.open("log.txt", FILE_WRITE);
  DateTime now;

digitalWrite(greenLEDpin, HIGH);
  // fetch the time
  now = RTC.now();
  logfile.print(now.year(), DEC);
  logfile.print("/");
  logfile.print(now.month(), DEC);
  logfile.print("/");
  logfile.print(now.day(), DEC);
  logfile.print(",");
  logfile.print(now.hour(), DEC);
  logfile.print(":");
  logfile.print(now.minute(), DEC);
  logfile.print(":");
  logfile.print(now.second(), DEC);
  logfile.print(",");
#if ECHO_TO_SERIAL
  Serial.print(now.year(), DEC);
  Serial.print("/");
  Serial.print(now.month(), DEC);
  Serial.print("/");
  Serial.print(now.day(), DEC);
  Serial.print(",");
  Serial.print(now.hour(), DEC);
  Serial.print(":");
  Serial.print(now.minute(), DEC);
  Serial.print(":");
  Serial.print(now.second(), DEC);
  Serial.print(",");
#endif //ECHO_TO_SERIAL

  // Reading temperature or humidity takes about 250 milliseconds!
  // Sensor readings may also be up to 2 seconds 'old' (its a very slow sensor)
  float h = dht.readHumidity();
  
  // Read temperature as Celsius
  float t = dht.readTemperature();
  
  // Check if any reads failed and exit early (to try again).
  //if (isnan(h) || isnan(t) || isnan(f)) {
    if (isnan(h) || isnan(t)) {
    Serial.println("Failed to read from DHT sensor!");
    return;
  }


  logfile.print(h); //log humidity
  logfile.print(",");
  logfile.println(t); //log temp (C)
   
#if ECHO_TO_SERIAL
  Serial.print(h); //display humidity
  Serial.print(",");
  Serial.println(t); //display temp (C)
#endif
  
  digitalWrite(greenLEDpin, LOW);
  
 logfile.close();

}


