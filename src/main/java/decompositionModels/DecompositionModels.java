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
package decompositionModels;

import org.jgrasstools.gears.libs.modules.JGTConstants;

import java.util.HashMap;
import java.util.Set;
import java.util.Map.Entry;

import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import decompositionModels.SimpleModelFactory;
import oms3.annotations.Author;
import oms3.annotations.Bibliography;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import oms3.annotations.Unit;

import org.geotools.feature.SchemaException;

import org.jgrasstools.gears.libs.modules.JGTModel;

// TODO: Auto-generated Javadoc
/**
 * The Class DecompositionModels.
 */
@Description("Cmputes the shortwave accounting for the cloudness")
@Documentation("")
@Author(name = "Giuseppe Formetta, Marialaura Bancheri", contact = "giuseppe.formetta@unitn.it")
@Keywords("Hydrology, Radiation, SkyviewFactor, Hillshade")
@Bibliography("Formetta (2013)")
@Label(JGTConstants.HYDROGEOMORPHOLOGY)
@Name("shortradbal")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
public class DecompositionModels extends JGTModel {

	@Description("Clearness index input Hashmap")
	@In
	@Unit("-")
	public HashMap<Integer, double[]> inClearnessIndexValues;

	@Description("The Hashmap with the time series of the SWRB meadured values")
	@In
	@Unit("W/m2")
	public HashMap<Integer, double[]> inSWRBMeasuredValues;

	@Description("The Hashmap with the time series of the direct shortwave radiation values")
	@In
	@Unit("W/m2")
	public HashMap<Integer, double[]> inSWRBDirectValues;

	@Description("The Hashmap with the time series of the diffuse shortwave radiation values")
	@In
	@Unit("W/m2")
	public HashMap<Integer, double[]> inSWRBDiffuseValues;

	@Description("String containing the name of the model: " + " Erbs; Reindl;Boland")
	@In
	public String model;

	Model DecMod;

	@Description("The output downwelling Hashmap")
	@Out
	public HashMap<Integer, double[]> outHMSWRBallSky = new HashMap<Integer, double[]>();;

	/**
	 * Process.
	 *
	 * @throws Exception the exception
	 */
	@Execute
	public void process() throws Exception {

		checkNull(inSWRBDirectValues);

		// reading the ID of all the stations
		Set<Entry<Integer, double[]>> entrySet = inSWRBDirectValues.entrySet();

		for (Entry<Integer, double[]> entry : entrySet) {

			Integer ID = entry.getKey();

			double SWRBMeasured = inSWRBMeasuredValues.get(ID)[0];

			double clearnessIndex = inClearnessIndexValues.get(ID)[0];

			if (isNovalue(clearnessIndex))
				clearnessIndex = 1;

			double SWRBdirect = inSWRBDirectValues.get(ID)[0];
			if (SWRBdirect < 0)
				SWRBdirect = 0;

			double SWRBdiffuse = inSWRBDiffuseValues.get(ID)[0];
			if (SWRBdiffuse < 0)
				SWRBdiffuse = 0;

			double SWRBallSky = 1. * SWRBdirect + 1. * SWRBdiffuse;

			if (!isNovalue(SWRBMeasured) & SWRBMeasured>0.) {

				double kd = computeKd(model, clearnessIndex);

				double cs = (SWRBdirect == 0) ? 0 : computeCs(kd, SWRBMeasured, SWRBdirect);

				double cd = (SWRBdiffuse == 0) ? 0 : computeCd(kd, SWRBMeasured, SWRBdiffuse);

				SWRBallSky = cd * SWRBdirect + cs * SWRBdiffuse;
			}

			/** Store results in Hashmaps */
			storeResult((Integer) ID, SWRBallSky);

		}

	}

	/**
	 * Compute kd which is the diffuse sky fraction coefficient (see also Helbig et
	 * al.2010).
	 *
	 * @param model          is the string containing the decomposition model
	 * @param clearnessIndex the clearness index value
	 * @return the double variable of the diffuse sky fraction coefficient
	 */
	private double computeKd(String model, double clearnessIndex) {
		DecMod = SimpleModelFactory.createModel(model, clearnessIndex);
		double kd = DecMod.kdValues();
		return kd = (kd > 1) ? 1 : kd;
	}

	/**
	 * Compute cs which is the correction coefficient for the direct radiation.
	 *
	 * @param kd           is the diffuse sky fraction coefficient
	 * @param SWRBMeasured is the measured shortwave
	 * @param SWRBdirect   is the direct shortwave
	 * @return the double value of the correction coefficient for the direct
	 *         radiation
	 */
	private double computeCs(double kd, double SWRBMeasured, double SWRBdirect) {
		return (1 - kd) * Math.min(1., SWRBMeasured / SWRBdirect);
	}

	/**
	 * Compute cd the correction coefficient for the diffuse radiation.
	 *
	 * @param kd           is the diffuse sky fraction coefficient
	 * @param SWRBMeasured is the measured shortwave
	 * @param SWRBdiffuse  is the diffuse shortwave
	 * @return the double value of the correction coefficient for the diffuse
	 *         radiation
	 */
	private double computeCd(double kd, double SWRBMeasured, double SWRBdiffuse) {

		return kd * Math.min(1., SWRBMeasured / SWRBdiffuse);
	}

	/**
	 * Store result.
	 *
	 * @param ID         the id of the station
	 * @param SWRBallSky the SWRB computed at all sky conditions
	 * @throws SchemaException the schema exception
	 */
	private void storeResult(int ID, double SWRBallSky) throws SchemaException {

		outHMSWRBallSky.put(ID, new double[] { SWRBallSky });
	}

}
