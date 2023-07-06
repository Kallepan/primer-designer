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

    // add title
    svg.append('text')
      .attr('x', width / 2)
      .attr('y', margins.TOP / 3)
      .attr('text-anchor', 'middle')
      .style('fill', 'white')
      .text(this.getTitle());

    // add primers
    svg.selectAll('primers')
      .data(this.regionData.primers)
      .enter()
      .append('rect')
      .attr('x', (primer) => x(primer.x1))
      .attr('y', (primer) => y(primer.pool))
      .attr('width', (primer) => x(primer.x2) - x(primer.x1))
      .attr('height', 10)
      .attr('transform', `translate(${margins.LEFT}, ${margins.TOP})`)
      .style('fill', 'steelblue');
    
    const div = d3.select('body').append('div')
      .attr('class', 'tooltip')
      .style('background-color', 'white')
      .style('font-size', '12px')
      .style('border', 'solid 1px #313639')
      .style('border-radius', '8px')
      .style('padding', '.25rem')
      .style('position', 'absolute')
      .style('text-align', 'center')
      .style('opacity', 0);

    // add amplicons
    svg.selectAll('amplicons')
      .data(this.regionData.amplicons)
      .enter()
      .append('rect')
      .attr('x', (amplicon) => x(amplicon.x1))
      .attr('y', (amplicon) => y(amplicon.pool))
      .attr('width', (amplicon) => x(amplicon.x2) - x(amplicon.x1))
      .attr('height', 10)
      .attr('transform', `translate(${margins.LEFT}, ${margins.TOP})`)
      .on('mouseover', (event, amplicon) => {
        d3.select(event.currentTarget).transition()
          .duration(50)
          .attr('opacity', '.85');
        div.transition()
          .duration(50)
          .style('opacity', 1);
        div.html(`Amplicon: ${amplicon.name}<br>Start: ${amplicon.x1}<br>End: ${amplicon.x2}`)
          .style('left', `${event.pageX + 10}px`)
          .style('top', `${event.pageY - 15}px`);
      })
      .on('mouseout', (event, amplicon) => {
        d3.select(event.currentTarget).transition()
          .duration(50)
          .attr('opacity', '1');
        div.transition()
          .duration(50)
          .style('opacity', 0);
      })
      .style('fill', 'red');

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