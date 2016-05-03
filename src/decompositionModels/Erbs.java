/*
 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
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


// TODO: Auto-generated Javadoc
/**
 * The Class Erbs implmentes Erbs [1982] decomposition model .
 */
public class Erbs implements Model{


	double cleranessIndex;
	
	double kd;


	/**
	 * Instantiates a new erbs model.
	 *
	 * @param cleranessIndex the cleraness index
	 */
	public Erbs(double cleranessIndex){

		this.cleranessIndex=cleranessIndex;
	}




	/* (non-Javadoc)
	 * @see decompositionModels.Model#kdValues()
	 */
	@Override
	public double kdValues() {
		
		if (cleranessIndex < 0.22) {
			kd = 1.0 - 0.09 * cleranessIndex;
		} else if (cleranessIndex < 0.80) {
			kd = 0.9511 - 0.1604 * cleranessIndex + 4.388 * Math.pow(cleranessIndex, 2.0) - 16.638
					* Math.pow(cleranessIndex, 3.0) + 12.336 * Math.pow(cleranessIndex, 4.0);
		} else {
			kd = 0.165;
		}
		return kd;
	}


}
