import { ComponentFixture, TestBed } from '@angular/core/testing';

import { SimplePlotComponent } from './simple-plot.component';

describe('PlotComponent', () => {
  let component: SimplePlotComponent;
  let fixture: ComponentFixture<SimplePlotComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [SimplePlotComponent]
    });
    fixture = TestBed.createComponent(SimplePlotComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
