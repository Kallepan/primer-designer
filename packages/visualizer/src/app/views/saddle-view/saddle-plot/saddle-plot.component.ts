import { AfterViewInit, Component, Input } from '@angular/core';
import * as d3 from 'd3';
import { SaddlePlotData } from '../types';
import { CONFIG } from 'src/app/config';

@Component({
  selector: 'app-saddle-plot',
  templateUrl: './saddle-plot.component.html',
  styleUrls: ['./saddle-plot.component.scss']
})
export class SaddlePlotComponent implements AfterViewInit {
  @Input() saddleData: SaddlePlotData = {} as SaddlePlotData;

  private _generatePlot(): void {
    const losses = this.saddleData.losses;
    const margins = CONFIG.SADDLE_VIEW.PLOT_MARGINS;
    const width = CONFIG.SADDLE_VIEW.PLOT_WIDTH;
    const height = CONFIG.SADDLE_VIEW.PLOT_HEIGHT;
    const id = "#" + this.getId();

    // D3 code here
    const svg = d3.select(id)
      .attr('width', width)
      .attr('height', height);

    // Compute scales and axes
    const xScale = d3.scaleLinear()
      .domain([1, losses.length])
      .range([margins.LEFT, width - margins.RIGHT]);

    const yScale = d3.scaleLinear()
      .domain([Math.min(...losses), Math.max(...losses)])
      .range([height - margins.BOTTOM, margins.TOP]);

    const xAxis = d3.axisBottom(xScale);
    const yAxis = d3.axisLeft(yScale);

    // Add axes
    svg.append('g')
      .attr('transform', `translate(0, ${height - margins.BOTTOM})`)
      .call(xAxis);

    svg.append('g')
      .attr('transform', `translate(${margins.LEFT}, 0)`)
      .call(yAxis);

    // Add axis labels
    svg.append('text')
      .attr('transform', `translate(${width / 2}, ${height})`)
      .style('text-anchor', 'middle')
      .style('fill', 'white')
      .text('Iteration');

    svg.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('y', -40)
      .attr('x', 0 - (height / 2))
      .attr('dy', '1em')
      .style('text-anchor', 'middle')
      .style('fill', 'white')
      .text('Loss');

    // Add line
    const line = d3.line<number>()
      .x((_, index) => xScale(index + 1))
      .y((loss) => yScale(loss));

    svg.append('path')
      .datum(losses)
      .attr('fill', 'none')
      .attr('stroke', '#b57efc')
      .attr('stroke-width', 1.5)
      .attr('d', line);
  }

  ngAfterViewInit(): void {
    this._generatePlot();
  }

  getId(): string {
    const name = this.saddleData.name;
    return `plot-${name}`;
  }
}
