import { AfterViewInit, Component, Input } from '@angular/core';
import { SimplifiedRegionData, SimplifiedPrimerData } from '../types';
import * as d3 from 'd3';
import { CONFIG } from 'src/app/config';

@Component({
  selector: 'app-plot',
  templateUrl: './plot.component.html',
  styleUrls: ['./plot.component.scss']
})
export class PlotComponent implements AfterViewInit {
  @Input() regionData: SimplifiedRegionData = {} as SimplifiedRegionData;

  private _generatePlot(): void {
    const pools = [...new Set(this.regionData.primers.map((primer) => primer.pool))];
    const margins = CONFIG.SIMPLIFIED_VIEW.PLOT_MARGINS;
    const width = CONFIG.SIMPLIFIED_VIEW.PLOT_WIDTH;
    const height = pools.length * CONFIG.SIMPLIFIED_VIEW.PRIMER_OFFSET + margins.TOP + margins.BOTTOM;
    const name = this.regionData.name;
    const id = `#${name}-plot`;
    
    // D3 code here
    const svg = d3.select(id)
      .attr('width', width)
      .attr('height', height);

    // x-axis
    const x = d3.scaleLinear()
      .domain([this.regionData.start, this.regionData.end])
      .range([margins.LEFT, width - margins.RIGHT]);
    
    // add x-axis
    svg.append('g')
      .attr('transform', `translate(0, ${height - margins.BOTTOM})`)
      .call(d3.axisBottom(x));
    
    // add x-axis label
    svg.append('text')
      .attr('transform', `translate(${width / 2}, ${height})`)
      .style('text-anchor', 'middle')
      .style('fill', 'white')
      .text('Position');

    // y-axis
    // Get unique pool names
    const range = pools.map((_, index) => index * CONFIG.SIMPLIFIED_VIEW.PRIMER_OFFSET);
    const y = d3.scaleOrdinal(pools, range);

    // add y-axis
    svg.append('g')
      .attr('transform', `translate(${margins.LEFT}, ${margins.TOP})`)
      .call(d3.axisLeft(y));
    
    // add y-axis label
    svg.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('y', 0)
      .attr('x', 0 - (height / 2))
      .attr('dy', '1em')
      .style('text-anchor', 'middle')
      .style('fill', 'white')
      .text('Pool');
  }

  getTitle(): string {
    return `${this.regionData?.name} (${this.regionData?.start}-${this.regionData?.end})`;
  }

  getId(): string {
    const name = this.regionData?.name;
    return `${name}-plot`;
  }

  ngAfterViewInit(): void {
    this._generatePlot();
  }
}
