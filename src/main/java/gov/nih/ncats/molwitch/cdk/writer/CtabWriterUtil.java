/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2020.
 *
 * This work is free software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License as published by the Free Software Foundation;
 * either version 2.1 of the License, or (at your option) any later version.
 *
 * This work is distributed in the hope that it will be useful, but without any warranty;
 * without even the implied warranty of merchantability or fitness for a particular purpose.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this library;
 *  if not, write to:
 *
 *  the Free Software Foundation, Inc.
 *  59 Temple Place, Suite 330
 *  Boston, MA 02111-1307 USA
 */

package gov.nih.ncats.molwitch.cdk.writer;

import gov.nih.ncats.molwitch.cdk.CdkUtil;
import gov.nih.ncats.molwitch.io.ChemFormat;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import java.util.Properties;
import java.util.function.Function;

public class CtabWriterUtil {
    public static Function<IAtomContainer, IAtomContainer> handleMolSpec(ChemFormat.MolFormatSpecification spec, Properties customSettings){
        ChemFormat.KekulizationEncoding kekulizationEncoding = spec.getKekulization();
//			 System.out.println("kekulizeEncoding = " + kekulizationEncoding);
        Function<IAtomContainer, IAtomContainer> adapter = c -> {
            IAtomContainer returnedContainer = c;
            boolean isAromatic = false;
            for (IBond bond : c.bonds()) {
                if (bond.isAromatic()) {
                    isAromatic = true;
                    break;
                }
            }
            switch (kekulizationEncoding) {

                case FORCE_AROMATIC: {
                    customSettings.setProperty("WriteAromaticBondTypes", Boolean.TRUE.toString());
                    //if we force aromaticity we need to set the aromatic flags

                    if (!isAromatic) {
                        try {
                            returnedContainer = returnedContainer.clone();
                            CdkUtil.aromatize(returnedContainer);
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }
                break;
                case KEKULE: {
                    customSettings.setProperty("WriteAromaticBondTypes", Boolean.FALSE.toString());
                    if(isAromatic){
                        try {
                            returnedContainer = CdkUtil.kekulizeIfNeeded(returnedContainer, true);
                        } catch (CDKException e) {
                            e.printStackTrace();
                        }
                    }
                }
                break;
                case AS_IS: {

                    customSettings.setProperty("WriteAromaticBondTypes", Boolean.toString(isAromatic));

                }
            }
            return returnedContainer;
        };


        Function<IAtomContainer, IAtomContainer> coordinateAdapter = Function.identity();
        switch (spec.getCoordinateOptions()) {
            case AS_IS: break;
            case FORCE_2D:
                customSettings.setProperty(
                        "OptForceWriteAs2DCoordinates", Boolean.toString(true)
                );

                coordinateAdapter = c -> {
                    //force all coords as 2d
                    boolean mustChange = false;
                    for (IAtom a : c.atoms()) {
                        if (a.getPoint2d() == null && a.getPoint3d() != null) {
                            mustChange = true;
                        }
                    }
                    if (mustChange) {
                        IAtomContainer clone;
                        try {
                            clone = c.clone();
                        } catch (CloneNotSupportedException e) {
                            return c;
                        }
                        for (IAtom a : clone.atoms()) {
                            Point3d points = a.getPoint3d();
                            if (points != null) {
                                a.setPoint2d(new Point2d(points.x, points.y));
                            }
                        }
                        return clone;
                    } else {
                        return c;
                    }
                };
                break;

            case FORCE_3D:
                customSettings.setProperty(
                        "OptForceWriteAs2DCoordinates", Boolean.toString(false)
                );
                coordinateAdapter = c -> {
                    //force all coords as 2d
                    boolean mustChange = false;
                    for (IAtom a : c.atoms()) {
                        if (a.getPoint3d() == null && a.getPoint2d() != null) {
                            mustChange = true;
                        }
                    }
                    if (mustChange) {
                        IAtomContainer clone;
                        try {
                            clone = c.clone();
                        } catch (CloneNotSupportedException e) {
                            return c;
                        }
                        for (IAtom a : clone.atoms()) {
                            Point2d points = a.getPoint2d();

                            if(points !=null){

                                a.setPoint3d(new Point3d(points.x, points.y, 0));
                            }
                        }
                        return clone;
                    } else {
                        return c;
                    }
                };

                break;
        }
        if (adapter == null) {
            adapter = coordinateAdapter;
        } else {
            adapter = adapter.andThen(coordinateAdapter);
        }
        return adapter;
    }
}
