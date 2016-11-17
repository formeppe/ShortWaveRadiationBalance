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

import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;


import decompositionModels.DecompositionModels;
import org.junit.Test;

/**
 * Test the {@link SWRB} module.
 * 
 * @author Marialaura Bancheri 
 */
public class TestDecompositionModels  {

	@Test
	public void Test() throws Exception {

		String startDate = "2007-10-17 00:00" ;
		String endDate = "2007-10-18 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";


		String inPathToCI ="resources/Input/CI.csv";
		String inPathToSWRBMeasured ="resources/Input/SWRBMeasured.csv";
		String inPathToSWRBDirect ="resources/Output/Direct.csv";
		String inPathToSWRBDiffuse ="resources/Output/Diffuse.csv";
		
		

		OmsTimeSeriesIteratorReader ciReader = getTimeseriesReader(inPathToCI, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader SWRBMeasured = getTimeseriesReader(inPathToSWRBMeasured, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader SWRBDirect = getTimeseriesReader(inPathToSWRBDirect, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader SWRBDiffuse = getTimeseriesReader(inPathToSWRBDiffuse, fId, startDate, endDate, timeStepMinutes);

		String pathToSWRBallSky ="resources/Output/SWRBallSky.csv";


		OmsTimeSeriesIteratorWriter writerSWRBallSky = new OmsTimeSeriesIteratorWriter();


		writerSWRBallSky .file = pathToSWRBallSky;
		writerSWRBallSky .tStart = startDate;
		writerSWRBallSky .tTimestep = timeStepMinutes;
		writerSWRBallSky .fileNovalue="-9999";


		DecompositionModels decMod = new DecompositionModels();

	
				

		while( ciReader.doProcess  ) {

			ciReader.nextRecord();	
			HashMap<Integer, double[]> id2ValueMap = ciReader.outData;
			decMod.inClearnessIndexValues= id2ValueMap;

			
			SWRBMeasured.nextRecord();
			id2ValueMap = SWRBMeasured.outData;
			decMod.inSWRBMeasuredValues = id2ValueMap;
			
			SWRBDirect.nextRecord();
			id2ValueMap = SWRBDirect.outData;
			decMod.inSWRBDirectValues = id2ValueMap;
			
			SWRBDiffuse.nextRecord();
			id2ValueMap = SWRBDiffuse.outData;
			decMod.inSWRBDiffuseValues = id2ValueMap;

			decMod.model="Reindl";
			decMod.process();

			HashMap<Integer, double[]> outHMdirect = decMod.outHMSWRBallSky;



			writerSWRBallSky.inData = outHMdirect ;
			writerSWRBallSky.writeNextLine();

			if (pathToSWRBallSky != null) {
				writerSWRBallSky .close();
			}



		}

		ciReader.close();
		SWRBMeasured.close();
		SWRBDirect.close();
		SWRBDiffuse.close();


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
