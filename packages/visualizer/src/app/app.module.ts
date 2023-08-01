import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';

import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { MaterialModule } from './material/material.module';
import { FormsModule, ReactiveFormsModule } from '@angular/forms';
import { NavigationComponent } from './navigation/navigation.component';
import { HomeComponent } from './home/home.component';
import { APP_BASE_HREF, HashLocationStrategy, LocationStrategy } from '@angular/common';
import { DetailViewComponent } from './views/detail-view/detail-view.component';
import { SimpleViewComponent } from './views/simple-view/simple-view.component';
import { SimplePlotComponent } from './views/simple-view/simple-plot/simple-plot.component';
import { SaddleViewComponent } from './views/saddle-view/saddle-view.component';
import { SaddlePlotComponent } from './views/saddle-view/saddle-plot/saddle-plot.component';
import { StatisticsComponent } from './statistics/statistics.component';

@NgModule({
  declarations: [
    AppComponent,
    HomeComponent,
    NavigationComponent,
    HomeComponent,
    DetailViewComponent,
    SimpleViewComponent,
    SimplePlotComponent,
    SaddleViewComponent,
    SaddlePlotComponent,
    StatisticsComponent
  ],
  imports: [
    BrowserModule,
    AppRoutingModule,
    MaterialModule,
    BrowserAnimationsModule,
    FormsModule,
    ReactiveFormsModule
  ],
  providers: [
    { provide: APP_BASE_HREF, useValue: '/' },
    { provide: LocationStrategy, useClass: HashLocationStrategy }
  ],
  bootstrap: [AppComponent]
})
export class AppModule { }
