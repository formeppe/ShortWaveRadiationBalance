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
 * A simple design factory for creating Model objects.
 */
public class SimpleModelFactory {

	/**
	 * Creates a new Model object.
	 *
	 * @param type: the string containing the name of the model
	 * @param clearnessIndex
	 * @return the model chosen
	 */
	public static Model createModel(String type,double clearnessIndex){
		Model model=null;


		if (type.equals("Erbs")){
			model=new Erbs(clearnessIndex);

		}else if (type.equals("Reindl")){
			model=new Reindl(clearnessIndex);

			/**Swinbank [1963]*/
		}else if (type.equals("Boland")){
			model=new Boland(clearnessIndex);}			


		return model;


	}
}
