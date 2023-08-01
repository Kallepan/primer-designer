import { ComponentFixture, TestBed } from '@angular/core/testing';

import { SaddlePlotComponent } from './saddle-plot.component';

describe('SaddlePlotComponent', () => {
  let component: SaddlePlotComponent;
  let fixture: ComponentFixture<SaddlePlotComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [SaddlePlotComponent]
    });
    fixture = TestBed.createComponent(SaddlePlotComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
