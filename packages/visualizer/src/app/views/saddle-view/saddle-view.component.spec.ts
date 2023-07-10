import { ComponentFixture, TestBed } from '@angular/core/testing';

import { SaddleViewComponent } from './saddle-view.component';

describe('SaddleViewComponent', () => {
  let component: SaddleViewComponent;
  let fixture: ComponentFixture<SaddleViewComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [SaddleViewComponent]
    });
    fixture = TestBed.createComponent(SaddleViewComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
