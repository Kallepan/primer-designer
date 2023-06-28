import { Component } from '@angular/core';
import resultsData from '../assets/results.json';

@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.scss']
})
export class AppComponent {
  title = 'visualizer';
  pools = resultsData.map(result => result.pool_id);
}
