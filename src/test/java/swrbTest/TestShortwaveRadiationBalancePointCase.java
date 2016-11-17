/*
 * GNU GPL v3 License
 *
 * Copyright 2016 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package swrbTest;

import java.net.URISyntaxException;
import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;

import swrbPointCase.ShortwaveRadiationBalancePointCase;

import org.junit.Test;

/**
 * Test the {@link SWRB} module.
 * 
 * @author Marialaura Bancheri 
 */
public class TestShortwaveRadiationBalancePointCase {

	@Test
	public void Test() throws Exception {

		String startDate = "2007-10-17 00:00" ;
		String endDate = "2007-10-18 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";


		String inPathToAirT ="resources/Input/temperature.csv";
		String inPathToHumidity ="resources/Input/humidity.csv";

		OmsTimeSeriesIteratorReader airTReader = getTimeseriesReader(inPathToAirT, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader humidityReader = getTimeseriesReader(inPathToHumidity, fId, startDate, endDate, timeStepMinutes);

		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/stations.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsRasterReader DEMreader = new OmsRasterReader();
		DEMreader.file = "resources/Input/DEM.asc";
		DEMreader.fileNovalue = -9999.0;
		DEMreader.geodataNovalue = Double.NaN;
		DEMreader.process();
		GridCoverage2D pit = DEMreader.outRaster;

		OmsRasterReader SKYreader = new OmsRasterReader();
		SKYreader.file = "resources/Input/skyview.asc";
		SKYreader.fileNovalue = -9999.0;
		SKYreader.geodataNovalue = Double.NaN;
		SKYreader.process();
		GridCoverage2D skyviewfactor = SKYreader.outRaster;

		String pathToDirect= "resources/Output/Direct.csv";
		String pathToDiffuse= "resources/Output/Diffuse.csv";
		String pathToTopATM= "resources/Output/TopATM.csv";

		OmsTimeSeriesIteratorWriter writerDirect = new OmsTimeSeriesIteratorWriter();
		OmsTimeSeriesIteratorWriter writerDiffuse = new OmsTimeSeriesIteratorWriter();
		OmsTimeSeriesIteratorWriter writerTopAtm = new OmsTimeSeriesIteratorWriter();


		writerDirect.file = pathToDirect;
		writerDirect.tStart = startDate;
		writerDirect.tTimestep = timeStepMinutes;
		writerDirect.fileNovalue="-9999";

		writerDiffuse.file = pathToDiffuse;
		writerDiffuse.tStart = startDate;
		writerDiffuse.tTimestep = timeStepMinutes;
		writerDiffuse.fileNovalue="-9999";

		writerTopAtm.file = pathToTopATM;
		writerTopAtm.tStart = startDate;
		writerTopAtm.tTimestep = timeStepMinutes;
		writerTopAtm.fileNovalue="-9999";



		ShortwaveRadiationBalancePointCase SWRBPoint = new ShortwaveRadiationBalancePointCase();
		SWRBPoint.inStations = stationsFC;
		SWRBPoint.inDem = pit;
		SWRBPoint.inSkyview = skyviewfactor;
		SWRBPoint.tStartDate = startDate;
		SWRBPoint.fStationsid = "cat";
		SWRBPoint.doHourly= true;
		SWRBPoint.pCmO3=0.6;
		SWRBPoint.pAlphag=0.9;
		SWRBPoint.pVisibility=80;
	
				

		while( airTReader.doProcess  ) {

			airTReader.nextRecord();	
			HashMap<Integer, double[]> id2ValueMap = airTReader.outData;
			SWRBPoint.inTemperatureValues= id2ValueMap;

			
			humidityReader.nextRecord();
			id2ValueMap = humidityReader.outData;
			SWRBPoint.inHumidityValues = id2ValueMap;


			SWRBPoint.process();

			HashMap<Integer, double[]> outHMdirect = SWRBPoint.outHMdirect;
			HashMap<Integer, double[]> outHMdiffuse = SWRBPoint.outHMdiffuse;
			HashMap<Integer, double[]> outHMtop = SWRBPoint.outHMtopatm;


			writerDirect.inData = outHMdirect ;
			writerDirect.writeNextLine();

			if (pathToDirect != null) {
				writerDirect .close();
			}

			writerDiffuse.inData = outHMdiffuse ;
			writerDiffuse.writeNextLine();

			if (pathToDiffuse != null) {
				writerDiffuse .close();
			}	

			writerTopAtm.inData =outHMtop;
			writerTopAtm.writeNextLine();

			if (pathToTopATM != null) {
				writerTopAtm.close();
			}


		}

		airTReader.close();
		humidityReader.close();


	}

	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}

}
