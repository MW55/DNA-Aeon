import constraint.mCGR;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYDotRenderer;
import org.jfree.data.xy.*;
import org.jfree.ui.ApplicationFrame;

import java.awt.*;

/**
 * ConstrainedKaos  23.06.20
 *
 * @author Hannah Franziska Löchel
 */
public class Plot {


    /**
     * Method for plotting
     *
     * @param kaosCons which should be plotted
     * @param size     of dots
     */
    public static void plot(mCGR kaosCons, int size) {
        XYSeries points = new XYSeries("Points");

        for (int key : kaosCons.getRows().keySet()) {

            for (int col : kaosCons.getRows().get(key)) {
                points.add(col, key);
            }

        }


        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(points);
        XYDotRenderer dot = new XYDotRenderer();

        dot.setSeriesPaint(0, Color.BLACK);
        dot.setDotHeight(size);
        dot.setDotWidth(size);

        NumberAxis xaxis = new NumberAxis();
        NumberAxis yaxis = new NumberAxis();

        XYPlot plot = new XYPlot(dataset, xaxis, yaxis, dot);
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);

        plot.setBackgroundPaint(Color.WHITE);
        org.jfree.chart.JFreeChart chart = new JFreeChart(plot);
        ApplicationFrame cgr = new ApplicationFrame("ConstrainedKaos");
        cgr.setSize(500, 500);

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setMouseWheelEnabled(true);

        chart.removeLegend();
        chartPanel.setPopupMenu(null);
        cgr.setContentPane(chartPanel);
        cgr.setVisible(true);

    }
}
