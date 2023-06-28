import { ComponentFixture, TestBed } from '@angular/core/testing';

import { DetailedViewComponent } from './detailed-view.component';

describe('DetailedViewComponent', () => {
  let component: DetailedViewComponent;
  let fixture: ComponentFixture<DetailedViewComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [DetailedViewComponent]
    });
    fixture = TestBed.createComponent(DetailedViewComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
