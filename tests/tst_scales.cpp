#include <QtTest>

#include "rscales.h"

class TestScales : public QObject
{
    Q_OBJECT

private slots:

    void defaultsAreUnity();
    void settersAndGetters();
    void copyAndAssignment();
    void frameVariablesUseCorrectDimensions();
    void beamSectionScalesWithLength();
};

void TestScales::defaultsAreUnity()
{
    RScales scales;
    QCOMPARE(scales.getMetre(), 1.0);
    QCOMPARE(scales.getKilogram(), 1.0);
    QCOMPARE(scales.getSecond(), 1.0);
    QCOMPARE(scales.getAmpere(), 1.0);
    QCOMPARE(scales.getKelvin(), 1.0);
    QCOMPARE(scales.getCandela(), 1.0);
    QCOMPARE(scales.getMole(), 1.0);
}

void TestScales::settersAndGetters()
{
    RScales scales;
    scales.setMetre(2.0);
    scales.setKilogram(3.0);
    scales.setSecond(4.0);
    scales.setKelvin(5.0);

    QCOMPARE(scales.getMetre(), 2.0);
    QCOMPARE(scales.getKilogram(), 3.0);
    QCOMPARE(scales.getSecond(), 4.0);
    QCOMPARE(scales.getKelvin(), 5.0);
}

void TestScales::copyAndAssignment()
{
    RScales scales;
    scales.setMetre(2.5);
    scales.setMole(6.0);

    RScales copy(scales);
    QCOMPARE(copy.getMetre(), 2.5);
    QCOMPARE(copy.getMole(), 6.0);

    RScales assigned;
    assigned = scales;
    QCOMPARE(assigned.getMetre(), 2.5);
    QCOMPARE(assigned.getMole(), 6.0);

    // Modifying the original must not affect the copy.
    scales.setMetre(9.0);
    QCOMPARE(copy.getMetre(), 2.5);
}

void TestScales::frameVariablesUseCorrectDimensions()
{
    RScales scales;
    scales.setMetre(0.5);
    scales.setKilogram(3.0);
    scales.setSecond(2.0);

    QCOMPARE(scales.findScaleFactor(R_VARIABLE_ROTATION),1.0);
    QCOMPARE(scales.findScaleFactor(R_VARIABLE_ROTATION_Z),1.0);
    QCOMPARE(scales.findScaleFactor(R_VARIABLE_MOMENT),0.1875);
    QCOMPARE(scales.findScaleFactor(R_VARIABLE_MOMENT_X),0.1875);
}

void TestScales::beamSectionScalesWithLength()
{
    RModel model;
    model.setNNodes(2);
    model.getNode(1).set(2.0,0.0,0.0);
    RElement element(R_ELEMENT_BEAM1);
    element.setNodeId(0,0);
    element.setNodeId(1,1);
    model.addElement(element,true,0);

    RLine &line = model.getLine(0);
    line.setCrossArea(0.02);
    line.setSecondMomentY(2.0e-5);
    line.setSecondMomentZ(3.0e-5);
    line.setTorsionalConstant(4.0e-5);

    RScales scales;
    scales.setMetre(0.5);
    scales.downscale(model);
    QCOMPARE(model.getNode(1).getX(),1.0);
    QCOMPARE(line.getCrossArea(),0.005);
    QCOMPARE(line.getSecondMomentY(),1.25e-6);
    QCOMPARE(line.getSecondMomentZ(),1.875e-6);
    QCOMPARE(line.getTorsionalConstant(),2.5e-6);

    scales.upscale(model);
    QCOMPARE(model.getNode(1).getX(),2.0);
    QCOMPARE(line.getCrossArea(),0.02);
    QCOMPARE(line.getSecondMomentY(),2.0e-5);
    QCOMPARE(line.getSecondMomentZ(),3.0e-5);
    QCOMPARE(line.getTorsionalConstant(),4.0e-5);
}

QTEST_APPLESS_MAIN(TestScales)

#include "tst_scales.moc"
